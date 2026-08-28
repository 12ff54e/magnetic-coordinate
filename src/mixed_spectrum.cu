#include "magnetic_coordinate/mixed_spectrum.hpp"

#include <cuda_runtime.h>
#include <cufft.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
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
    explicit DeviceBuffer(std::size_t count) : count_(count) {
        if (count_ != 0) {
            check_cuda(cudaMalloc(&data_, count_ * sizeof(T)), "cudaMalloc");
        }
    }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    ~DeviceBuffer() { cudaFree(data_); }
    T* data() { return data_; }

   private:
    T* data_ = nullptr;
    std::size_t count_ = 0;
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

void validate_input(const BoozerMixedGrid& grid, int mmax, int nmax) {
    if (grid.source_ns <= grid.first_surface || grid.first_surface < 1 ||
        grid.ntheta < 2 || grid.nzeta < 1 || mmax < 0 || nmax < 0 ||
        2 * mmax >= grid.ntheta || 2 * nmax >= grid.nzeta) {
        throw std::invalid_argument("invalid mixed-grid spectrum dimensions");
    }
    const std::size_t count =
        static_cast<std::size_t>(grid.source_ns - grid.first_surface) *
        static_cast<std::size_t>(grid.ntheta) * grid.nzeta;
    if (grid.r.size() != count || grid.z.size() != count ||
        grid.nu.size() != count) {
        throw std::invalid_argument("mixed-grid spectrum extent mismatch");
    }
    for (std::size_t index = 0; index < count; ++index) {
        if (!std::isfinite(grid.r[index]) || !std::isfinite(grid.z[index]) ||
            !std::isfinite(grid.nu[index])) {
            throw std::invalid_argument(
                "mixed-grid spectrum input contains non-finite values");
        }
    }
}

}  // namespace

MixedGridSpectrum analyze_mixed_grid_gpu(const BoozerMixedGrid& grid,
                                         int mmax,
                                         int nmax) {
    validate_input(grid, mmax, nmax);
    const int surfaces = grid.source_ns - grid.first_surface;
    const int points = grid.ntheta * grid.nzeta;
    const std::size_t field_count =
        static_cast<std::size_t>(surfaces) * static_cast<std::size_t>(points);
    constexpr int fields = 3;
    std::vector<cufftDoubleComplex> host_values(
        static_cast<std::size_t>(fields) * field_count);
    for (std::size_t index = 0; index < field_count; ++index) {
        host_values[index] = {grid.r[index], 0.0};
        host_values[field_count + index] = {grid.z[index], 0.0};
        host_values[2 * field_count + index] = {grid.nu[index], 0.0};
    }

    DeviceBuffer<cufftDoubleComplex> device_values(host_values.size());
    check_cuda(cudaMemcpy(device_values.data(), host_values.data(),
                          host_values.size() * sizeof(cufftDoubleComplex),
                          cudaMemcpyHostToDevice),
               "copy mixed-grid fields to GPU");
    CufftPlan plan(grid.nzeta, grid.ntheta, fields * surfaces);
    check_cufft(cufftExecZ2Z(plan.get(), device_values.data(),
                             device_values.data(), CUFFT_FORWARD),
                "forward mixed-grid FFT");
    check_cuda(cudaMemcpy(host_values.data(), device_values.data(),
                          host_values.size() * sizeof(cufftDoubleComplex),
                          cudaMemcpyDeviceToHost),
               "copy mixed-grid spectra from GPU");

    MixedGridSpectrum result;
    result.source_ns = grid.source_ns;
    result.first_surface = grid.first_surface;
    result.mmax = mmax;
    result.nmax = nmax;
    const int mode_count = (2 * mmax + 1) * (2 * nmax + 1);
    result.m.reserve(static_cast<std::size_t>(mode_count));
    result.n.reserve(static_cast<std::size_t>(mode_count));
    for (int n = -nmax; n <= nmax; ++n) {
        for (int m = -mmax; m <= mmax; ++m) {
            result.m.push_back(m);
            result.n.push_back(n);
        }
    }
    const std::size_t coefficient_count =
        static_cast<std::size_t>(surfaces) *
        static_cast<std::size_t>(mode_count);
    result.r.resize(coefficient_count);
    result.z.resize(coefficient_count);
    result.nu.resize(coefficient_count);
    const double weight = 4.0 * std::numbers::pi * std::numbers::pi /
                          static_cast<double>(points);
    for (int surface = 0; surface < surfaces; ++surface) {
        for (int mode = 0; mode < mode_count; ++mode) {
            const int m = result.m[static_cast<std::size_t>(mode)];
            const int n = result.n[static_cast<std::size_t>(mode)];
            const int theta_index = m >= 0 ? m : grid.ntheta + m;
            const int zeta_index = n >= 0 ? n : grid.nzeta + n;
            const std::size_t source =
                static_cast<std::size_t>(surface) * points +
                static_cast<std::size_t>(zeta_index * grid.ntheta +
                                         theta_index);
            const std::size_t output =
                static_cast<std::size_t>(mode) +
                static_cast<std::size_t>(mode_count * surface);
            const auto convert = [weight](cufftDoubleComplex value) {
                return std::complex<double>(weight * value.x,
                                            weight * value.y);
            };
            result.r[output] = convert(host_values[source]);
            result.z[output] = convert(host_values[field_count + source]);
            result.nu[output] = convert(host_values[2 * field_count + source]);
        }
    }
    return result;
}

}  // namespace magnetic_coordinate
