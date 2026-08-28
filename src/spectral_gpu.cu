#include "magnetic_coordinate/spectral.hpp"

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace magnetic_coordinate {
namespace {

void check_cuda(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(status));
    }
}

void check_cufft(cufftResult status, const char* operation) {
    if (status != CUFFT_SUCCESS) {
        throw std::runtime_error(std::string(operation) +
                                 " failed with cuFFT status " +
                                 std::to_string(static_cast<int>(status)));
    }
}

template <typename T>
class DeviceBuffer {
   public:
    explicit DeviceBuffer(std::size_t count) {
        if (count != 0) {
            check_cuda(cudaMalloc(&data_, count * sizeof(T)), "cudaMalloc");
        }
    }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    ~DeviceBuffer() { cudaFree(data_); }
    T* data() { return data_; }

   private:
    T* data_ = nullptr;
};

class CufftPlan {
   public:
    CufftPlan(int nzeta, int ntheta, int batch) {
        int dimensions[2] = {nzeta, ntheta};
        check_cufft(cufftPlanMany(&plan_, 2, dimensions, nullptr, 1,
                                  nzeta * ntheta, nullptr, 1,
                                  nzeta * ntheta, CUFFT_Z2Z, batch),
                    "cufftPlanMany");
    }
    CufftPlan(const CufftPlan&) = delete;
    CufftPlan& operator=(const CufftPlan&) = delete;
    ~CufftPlan() {
        if (plan_ != 0) cufftDestroy(plan_);
    }
    cufftHandle get() const { return plan_; }

   private:
    cufftHandle plan_ = 0;
};

__device__ void add_complex(cufftDoubleComplex* destination,
                            double real,
                            double imaginary) {
    atomicAdd(&destination->x, real);
    atomicAdd(&destination->y, imaginary);
}

// Expand each folded (m>=0,n>=0) product-basis coefficient into its four
// signed complex modes. Atomic accumulation deliberately handles coincident
// zero and Nyquist modes without separate normalization cases.
__global__ void assemble_geometry_spectra_kernel(
    const double* __restrict__ families,
    const double* __restrict__ lambda_multiplier,
    cufftDoubleComplex* __restrict__ spectra,
    int ns,
    int mpol,
    int ntor,
    int nfp,
    int ntheta,
    int nzeta,
    int points,
    int family_size,
    int coefficients) {
    const int work = blockIdx.x * blockDim.x + threadIdx.x;
    if (work >= 4 * coefficients) return;
    const int coefficient = work / 4;
    const int quadrant = work - 4 * coefficient;
    const int sign_m = (quadrant & 1) == 0 ? 1 : -1;
    const int sign_n = (quadrant & 2) == 0 ? 1 : -1;
    const int surface = coefficient % ns;
    const int mode = coefficient / ns;
    const int m = mode / (ntor + 1);
    const int n = mode - m * (ntor + 1);
    if (m >= mpol) return;
    const int signed_m = sign_m * m;
    const int signed_n = sign_n * n;
    const int theta_index = signed_m >= 0 ? signed_m : ntheta + signed_m;
    const int zeta_index = signed_n >= 0 ? signed_n : nzeta + signed_n;
    const int spectrum_point = surface * points +
                               zeta_index * ntheta + theta_index;
    const int field_count = ns * points;
    constexpr double quarter = 0.25;

    const double rcc = families[coefficient];
    const double zsc = families[family_size + coefficient];
    const double lsc = families[2 * family_size + coefficient];
    const double rss = families[3 * family_size + coefficient];
    const double zcs = families[4 * family_size + coefficient];
    const double lcs = families[5 * family_size + coefficient];
    const double multiplier = lambda_multiplier[surface];
    const double r_real =
        quarter * (rcc - static_cast<double>(sign_m * sign_n) * rss);
    const double z_imaginary =
        -quarter * (static_cast<double>(sign_m) * zsc +
                    static_cast<double>(sign_n) * zcs);
    const double lambda_imaginary =
        -quarter * multiplier *
        (static_cast<double>(sign_m) * lsc +
         static_cast<double>(sign_n) * lcs);

    add_complex(&spectra[spectrum_point], r_real, 0.0);
    add_complex(&spectra[field_count + spectrum_point], 0.0, z_imaginary);
    add_complex(&spectra[2 * field_count + spectrum_point], 0.0,
                lambda_imaginary);
    // i*k times an imaginary lambda coefficient is purely real.
    add_complex(&spectra[3 * field_count + spectrum_point],
                -static_cast<double>(signed_m) * lambda_imaginary, 0.0);
    add_complex(&spectra[4 * field_count + spectrum_point],
                -static_cast<double>(signed_n * nfp) * lambda_imaginary, 0.0);
}

void validate_input(const CumesEquilibriumView& equilibrium,
                    const FluxNormalization& normalization) {
    if (equilibrium.ns < 3 || equilibrium.mpol < 1 || equilibrium.ntor < 0 ||
        equilibrium.nfp < 1 || equilibrium.ntheta < 1 ||
        equilibrium.nzeta < 1 ||
        equilibrium.mnmax != equilibrium.mpol * (equilibrium.ntor + 1) ||
        equilibrium.mpol - 1 > equilibrium.ntheta / 2 ||
        equilibrium.ntor > equilibrium.nzeta / 2 ||
        normalization.lambda_multiplier.size() !=
            static_cast<std::size_t>(equilibrium.ns)) {
        throw std::invalid_argument("invalid GPU spectral dimensions");
    }
    const std::size_t family_size =
        static_cast<std::size_t>(equilibrium.ns) * equilibrium.mnmax;
    for (const auto& family : equilibrium.families) {
        if (family.size() != family_size) {
            throw std::invalid_argument("GPU spectral family extent mismatch");
        }
        for (double value : family) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "GPU spectral family contains a non-finite value");
            }
        }
    }
    for (double value : normalization.lambda_multiplier) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "GPU lambda normalization contains a non-finite value");
        }
    }
}

}  // namespace

NativeAngularGeometry synthesize_native_geometry_gpu(
    const CumesEquilibriumView& equilibrium,
    const FluxNormalization& normalization) {
    validate_input(equilibrium, normalization);
    const std::size_t family_size =
        static_cast<std::size_t>(equilibrium.ns) * equilibrium.mnmax;
    std::vector<double> host_families(
        CumesEquilibrium::SPECTRAL_FAMILY_COUNT * family_size);
    for (std::size_t family = 0;
         family < CumesEquilibrium::SPECTRAL_FAMILY_COUNT; ++family) {
        std::copy(equilibrium.families[family].begin(),
                  equilibrium.families[family].end(),
                  host_families.begin() +
                      static_cast<std::ptrdiff_t>(family * family_size));
    }
    const int points = equilibrium.ntheta * equilibrium.nzeta;
    const int field_count = equilibrium.ns * points;
    constexpr int output_fields = 5;
    const std::size_t spectrum_count =
        static_cast<std::size_t>(output_fields) * field_count;
    DeviceBuffer<double> device_families(host_families.size());
    DeviceBuffer<double> device_multiplier(
        normalization.lambda_multiplier.size());
    DeviceBuffer<cufftDoubleComplex> device_spectra(spectrum_count);
    check_cuda(cudaMemcpy(device_families.data(), host_families.data(),
                          host_families.size() * sizeof(double),
                          cudaMemcpyHostToDevice),
               "copy spectral families to GPU");
    check_cuda(cudaMemcpy(device_multiplier.data(),
                          normalization.lambda_multiplier.data(),
                          normalization.lambda_multiplier.size() *
                              sizeof(double),
                          cudaMemcpyHostToDevice),
               "copy lambda normalization to GPU");
    check_cuda(cudaMemset(device_spectra.data(), 0,
                          spectrum_count * sizeof(cufftDoubleComplex)),
               "clear native geometry spectra");
    const int coefficients = static_cast<int>(family_size);
    const int work = 4 * coefficients;
    assemble_geometry_spectra_kernel<<<(work + 255) / 256, 256>>>(
        device_families.data(), device_multiplier.data(),
        device_spectra.data(), equilibrium.ns, equilibrium.mpol,
        equilibrium.ntor, equilibrium.nfp, equilibrium.ntheta,
        equilibrium.nzeta, points, coefficients, coefficients);
    check_cuda(cudaGetLastError(), "assemble native geometry spectra");
    CufftPlan plan(equilibrium.nzeta, equilibrium.ntheta,
                   output_fields * equilibrium.ns);
    check_cufft(cufftExecZ2Z(plan.get(), device_spectra.data(),
                             device_spectra.data(), CUFFT_INVERSE),
                "inverse native geometry FFT");

    std::vector<cufftDoubleComplex> host_output(spectrum_count);
    check_cuda(cudaMemcpy(host_output.data(), device_spectra.data(),
                          host_output.size() * sizeof(cufftDoubleComplex),
                          cudaMemcpyDeviceToHost),
               "copy native geometry from GPU");
    NativeAngularGeometry result;
    result.ns = equilibrium.ns;
    result.ntheta = equilibrium.ntheta;
    result.nzeta = equilibrium.nzeta;
    std::vector<double>* outputs[output_fields] = {
        &result.r, &result.z, &result.lambda, &result.lambda_theta,
        &result.lambda_zeta};
    for (int field = 0; field < output_fields; ++field) {
        outputs[field]->resize(static_cast<std::size_t>(field_count));
        for (int point = 0; point < field_count; ++point) {
            const auto value = host_output[static_cast<std::size_t>(
                field * field_count + point)];
            if (!std::isfinite(value.x) || !std::isfinite(value.y) ||
                std::abs(value.y) >
                    1.0e-11 * std::max(1.0, std::abs(value.x))) {
                throw std::runtime_error(
                    "inverse geometry FFT did not produce a real field");
            }
            (*outputs[field])[static_cast<std::size_t>(point)] = value.x;
        }
    }
    return result;
}

}  // namespace magnetic_coordinate
