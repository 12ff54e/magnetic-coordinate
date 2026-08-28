#include "magnetic_coordinate/boozer_fft.hpp"

#include <cuda_runtime.h>
#include <cufft.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace magnetic_coordinate {
namespace {

constexpr double PI = 3.141592653589793238462643383279502884;

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
    explicit DeviceBuffer(std::size_t count) : count_(count) {
        if (count_ != 0) {
            check_cuda(cudaMalloc(&data_, count_ * sizeof(T)), "cudaMalloc");
        }
    }

    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    DeviceBuffer(DeviceBuffer&& other) noexcept
        : data_(std::exchange(other.data_, nullptr)),
          count_(std::exchange(other.count_, 0)) {}
    DeviceBuffer& operator=(DeviceBuffer&&) = delete;

    ~DeviceBuffer() { cudaFree(data_); }

    T* data() { return data_; }
    const T* data() const { return data_; }
    std::size_t size() const { return count_; }

   private:
    T* data_ = nullptr;
    std::size_t count_ = 0;
};

class CufftPlan {
   public:
    CufftPlan(int nzeta, int ntheta, int batch) {
        int dimensions[2] = {nzeta, ntheta};
        check_cufft(
            cufftPlanMany(&plan_, 2, dimensions, nullptr, 1, nzeta * ntheta,
                          nullptr, 1, nzeta * ntheta, CUFFT_Z2Z, batch),
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

__global__ void capture_zero_modes_kernel(
    const cufftDoubleComplex* __restrict__ spectrum,
    double* __restrict__ raw_zero,
    double* __restrict__ integral_zero,
    int points,
    int surfaces,
    int* __restrict__ status) {
    const int surface = blockIdx.x * blockDim.x + threadIdx.x;
    if (surface >= surfaces) return;
    const cufftDoubleComplex zero = spectrum[surface * points];
    if (!isfinite(zero.x) || !isfinite(zero.y) || zero.x == 0.0) {
        atomicCAS(status, 0, 1);
        return;
    }
    raw_zero[surface] = zero.x;
    integral_zero[surface] =
        (4.0 * PI * PI / static_cast<double>(points)) * zero.x;
}

__global__ void solve_shift_modes_kernel(
    cufftDoubleComplex* __restrict__ spectrum,
    const double* __restrict__ raw_zero,
    const double* __restrict__ iota,
    int points,
    int surfaces,
    int ntheta,
    int nzeta,
    int nfp,
    double tolerance,
    int* __restrict__ status) {
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    const int total = surfaces * points;
    if (index >= total) return;
    const int surface = index / points;
    const int local = index - surface * points;
    const int theta_index = local % ntheta;
    const int zeta_index = local / ntheta;
    if (local == 0) {
        spectrum[index] = cufftDoubleComplex{0.0, 0.0};
        return;
    }

    const bool theta_nyquist = ntheta % 2 == 0 && theta_index == ntheta / 2;
    const bool zeta_nyquist = nzeta % 2 == 0 && zeta_index == nzeta / 2;
    if (theta_nyquist || zeta_nyquist) {
        spectrum[index] = cufftDoubleComplex{0.0, 0.0};
        return;
    }

    const int m =
        theta_index <= ntheta / 2 ? theta_index : theta_index - ntheta;
    const int n = zeta_index <= nzeta / 2 ? zeta_index : zeta_index - nzeta;
    const double denominator =
        iota[surface] * static_cast<double>(m) + static_cast<double>(n * nfp);
    const double denominator_scale =
        fabs(iota[surface]) * static_cast<double>(ntheta / 2) +
        static_cast<double>(nfp * (nzeta / 2));
    const cufftDoubleComplex value = spectrum[index];
    if (fabs(denominator) <= tolerance * fmax(1.0, denominator_scale)) {
        if (hypot(value.x, value.y) > tolerance * fabs(raw_zero[surface])) {
            atomicCAS(status, 0, 2);
        }
        spectrum[index] = cufftDoubleComplex{0.0, 0.0};
        return;
    }

    const double scale = raw_zero[surface] * denominator;
    spectrum[index] = cufftDoubleComplex{value.y / scale, -value.x / scale};
}

void validate_input(const PestGrid& pest,
                    std::span<const double> iota,
                    int nfp,
                    double tolerance) {
    if (pest.source_ns < 2 || pest.first_surface < 1 || pest.ntheta < 2 ||
        pest.nzeta < 1 || nfp < 1 || !(tolerance > 0.0) ||
        !std::isfinite(tolerance)) {
        throw std::invalid_argument("invalid Boozer FFT configuration");
    }
    const int surfaces = pest.source_ns - pest.first_surface;
    const std::size_t points =
        static_cast<std::size_t>(pest.ntheta) * pest.nzeta;
    const std::size_t count = static_cast<std::size_t>(surfaces) * points;
    if (pest.b2.size() != count || pest.sqrtg.size() != count ||
        iota.size() != static_cast<std::size_t>(pest.source_ns)) {
        throw std::invalid_argument("Boozer FFT input extent mismatch");
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (!(pest.b2[index] > 0.0) || !std::isfinite(pest.b2[index]) ||
            !std::isfinite(pest.sqrtg[index])) {
            throw std::invalid_argument(
                "Boozer FFT input contains a nonphysical field");
        }
    }
    for (int surface = pest.first_surface; surface < pest.source_ns;
         ++surface) {
        if (!std::isfinite(iota[static_cast<std::size_t>(surface)])) {
            throw std::invalid_argument("Boozer FFT iota is non-finite");
        }
    }
}

}  // namespace

BoozerShift solve_boozer_shift_gpu(const PestGrid& pest,
                                   std::span<const double> iota,
                                   int nfp,
                                   double resonance_tolerance) {
    validate_input(pest, iota, nfp, resonance_tolerance);
    const int surfaces = pest.source_ns - pest.first_surface;
    const int points = pest.ntheta * pest.nzeta;
    const std::size_t count =
        static_cast<std::size_t>(surfaces) * static_cast<std::size_t>(points);

    std::vector<cufftDoubleComplex> host_spectrum(count);
    for (std::size_t index = 0; index < count; ++index) {
        host_spectrum[index] =
            cufftDoubleComplex{pest.b2[index] * pest.sqrtg[index], 0.0};
    }
    std::vector<double> surface_iota(static_cast<std::size_t>(surfaces));
    for (int surface = 0; surface < surfaces; ++surface) {
        surface_iota[static_cast<std::size_t>(surface)] =
            iota[static_cast<std::size_t>(surface + pest.first_surface)];
    }

    DeviceBuffer<cufftDoubleComplex> device_spectrum(count);
    DeviceBuffer<double> device_iota(static_cast<std::size_t>(surfaces));
    DeviceBuffer<double> device_raw_zero(static_cast<std::size_t>(surfaces));
    DeviceBuffer<double> device_integral_zero(
        static_cast<std::size_t>(surfaces));
    DeviceBuffer<int> device_status(1);
    check_cuda(
        cudaMemcpy(device_spectrum.data(), host_spectrum.data(),
                   count * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice),
        "copy q to GPU");
    check_cuda(cudaMemcpy(device_iota.data(), surface_iota.data(),
                          surface_iota.size() * sizeof(double),
                          cudaMemcpyHostToDevice),
               "copy iota to GPU");
    check_cuda(cudaMemset(device_status.data(), 0, sizeof(int)),
               "clear Boozer FFT status");

    CufftPlan plan(pest.nzeta, pest.ntheta, surfaces);
    check_cufft(cufftExecZ2Z(plan.get(), device_spectrum.data(),
                             device_spectrum.data(), CUFFT_FORWARD),
                "forward Boozer FFT");
    capture_zero_modes_kernel<<<(surfaces + 127) / 128, 128>>>(
        device_spectrum.data(), device_raw_zero.data(),
        device_integral_zero.data(), points, surfaces, device_status.data());
    check_cuda(cudaGetLastError(), "capture b2j00");

    const int blocks = static_cast<int>((count + 255) / 256);
    solve_shift_modes_kernel<<<blocks, 256>>>(
        device_spectrum.data(), device_raw_zero.data(), device_iota.data(),
        points, surfaces, pest.ntheta, pest.nzeta, nfp, resonance_tolerance,
        device_status.data());
    check_cuda(cudaGetLastError(), "solve Boozer shift modes");
    check_cufft(cufftExecZ2Z(plan.get(), device_spectrum.data(),
                             device_spectrum.data(), CUFFT_INVERSE),
                "inverse Boozer FFT");

    int status = 0;
    check_cuda(cudaMemcpy(&status, device_status.data(), sizeof(int),
                          cudaMemcpyDeviceToHost),
               "copy Boozer FFT status");
    if (status == 1) {
        throw std::domain_error("B^2 sqrt(g_p) has a zero or invalid mode 00");
    }
    if (status == 2) {
        throw std::domain_error(
            "B^2 sqrt(g_p) has an unresolved resonant Fourier mode");
    }

    check_cuda(
        cudaMemcpy(host_spectrum.data(), device_spectrum.data(),
                   count * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost),
        "copy Boozer shift from GPU");
    BoozerShift result;
    result.source_ns = pest.source_ns;
    result.first_surface = pest.first_surface;
    result.ntheta = pest.ntheta;
    result.nzeta = pest.nzeta;
    result.nu.resize(count);
    result.b2j00.resize(static_cast<std::size_t>(surfaces));
    check_cuda(cudaMemcpy(result.b2j00.data(), device_integral_zero.data(),
                          result.b2j00.size() * sizeof(double),
                          cudaMemcpyDeviceToHost),
               "copy b2j00 from GPU");
    for (std::size_t index = 0; index < count; ++index) {
        const double real = host_spectrum[index].x;
        const double imaginary = host_spectrum[index].y;
        if (!std::isfinite(real) || !std::isfinite(imaginary) ||
            std::abs(imaginary) > 1.0e-10 * std::max(1.0, std::abs(real))) {
            throw std::runtime_error(
                "inverse Boozer FFT did not produce a finite real field");
        }
        result.nu[index] = real;
    }
    return result;
}

}  // namespace magnetic_coordinate
