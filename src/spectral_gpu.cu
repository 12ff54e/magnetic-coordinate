#include "magnetic_coordinate/spectral.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace magnetic_coordinate {
namespace {

constexpr double TWO_PI = 6.283185307179586476925286766559005768;

void check_cuda(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string(operation) + ": " +
                                 cudaGetErrorString(status));
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

__global__ void synthesize_geometry_kernel(
    const double* __restrict__ families,
    const double* __restrict__ lambda_multiplier,
    double* __restrict__ output,
    int ns,
    int mpol,
    int ntor,
    int nfp,
    int ntheta,
    int nzeta,
    int points,
    int family_size,
    int total) {
    const int point = blockIdx.x * blockDim.x + threadIdx.x;
    if (point >= total) return;
    const int surface = point / points;
    const int angular = point - surface * points;
    const int theta_index = angular % ntheta;
    const int zeta_index = angular / ntheta;
    const double theta = TWO_PI * static_cast<double>(theta_index) /
                         static_cast<double>(ntheta);
    const double zeta = TWO_PI * static_cast<double>(zeta_index) /
                        static_cast<double>(nzeta);
    const double multiplier = lambda_multiplier[surface];
    double r = 0.0;
    double z = 0.0;
    double lambda = 0.0;
    double lambda_theta = 0.0;
    double lambda_zeta = 0.0;
    for (int m = 0; m < mpol; ++m) {
        double sin_m = 0.0;
        double cos_m = 0.0;
        sincos(static_cast<double>(m) * theta, &sin_m, &cos_m);
        for (int n = 0; n <= ntor; ++n) {
            double sin_n = 0.0;
            double cos_n = 0.0;
            sincos(static_cast<double>(n) * zeta, &sin_n, &cos_n);
            const int mode = m * (ntor + 1) + n;
            const int coefficient = surface + mode * ns;
            const double rcc = families[coefficient];
            const double zsc = families[family_size + coefficient];
            const double lsc = families[2 * family_size + coefficient];
            const double rss = families[3 * family_size + coefficient];
            const double zcs = families[4 * family_size + coefficient];
            const double lcs = families[5 * family_size + coefficient];
            r += rcc * cos_m * cos_n + rss * sin_m * sin_n;
            z += zsc * sin_m * cos_n + zcs * cos_m * sin_n;
            lambda += multiplier *
                      (lsc * sin_m * cos_n + lcs * cos_m * sin_n);
            lambda_theta +=
                multiplier * static_cast<double>(m) *
                (lsc * cos_m * cos_n - lcs * sin_m * sin_n);
            lambda_zeta +=
                multiplier * static_cast<double>(n * nfp) *
                (-lsc * sin_m * sin_n + lcs * cos_m * cos_n);
        }
    }
    output[point] = r;
    output[total + point] = z;
    output[2 * total + point] = lambda;
    output[3 * total + point] = lambda_theta;
    output[4 * total + point] = lambda_zeta;
}

void validate_input(const CumesEquilibrium& equilibrium,
                    const FluxNormalization& normalization) {
    if (equilibrium.ns < 3 || equilibrium.mpol < 1 || equilibrium.ntor < 0 ||
        equilibrium.nfp < 1 || equilibrium.ntheta < 1 ||
        equilibrium.nzeta < 1 ||
        equilibrium.mnmax != equilibrium.mpol * (equilibrium.ntor + 1) ||
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
    const CumesEquilibrium& equilibrium,
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
    const int total = equilibrium.ns * points;
    DeviceBuffer<double> device_families(host_families.size());
    DeviceBuffer<double> device_multiplier(
        normalization.lambda_multiplier.size());
    DeviceBuffer<double> device_output(static_cast<std::size_t>(5 * total));
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
    synthesize_geometry_kernel<<<(total + 255) / 256, 256>>>(
        device_families.data(), device_multiplier.data(),
        device_output.data(), equilibrium.ns, equilibrium.mpol,
        equilibrium.ntor, equilibrium.nfp, equilibrium.ntheta,
        equilibrium.nzeta, points, static_cast<int>(family_size), total);
    check_cuda(cudaGetLastError(), "synthesize native geometry");
    std::vector<double> host_output(static_cast<std::size_t>(5 * total));
    check_cuda(cudaMemcpy(host_output.data(), device_output.data(),
                          host_output.size() * sizeof(double),
                          cudaMemcpyDeviceToHost),
               "copy native geometry from GPU");

    NativeAngularGeometry result;
    result.ns = equilibrium.ns;
    result.ntheta = equilibrium.ntheta;
    result.nzeta = equilibrium.nzeta;
    const auto copy_field = [&](int field) {
        const auto first = host_output.begin() +
                           static_cast<std::ptrdiff_t>(field * total);
        return std::vector<double>(first, first + total);
    };
    result.r = copy_field(0);
    result.z = copy_field(1);
    result.lambda = copy_field(2);
    result.lambda_theta = copy_field(3);
    result.lambda_zeta = copy_field(4);
    return result;
}

}  // namespace magnetic_coordinate
