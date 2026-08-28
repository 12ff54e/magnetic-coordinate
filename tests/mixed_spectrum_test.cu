#include "magnetic_coordinate/mixed_spectrum.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

namespace {

void expect_near(std::complex<double> actual,
                 std::complex<double> expected,
                 double tolerance,
                 const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

std::size_t find_mode(const magnetic_coordinate::MixedGridSpectrum& spectrum,
                      int wanted_m,
                      int wanted_n) {
    for (std::size_t mode = 0; mode < spectrum.m.size(); ++mode) {
        if (spectrum.m[mode] == wanted_m && spectrum.n[mode] == wanted_n) {
            return mode;
        }
    }
    throw std::runtime_error("requested test mode is absent");
}

void test_integral_normalization() {
    constexpr int ntheta = 16;
    constexpr int nzeta = 12;
    magnetic_coordinate::BoozerMixedGrid grid;
    grid.source_ns = 2;
    grid.ntheta = ntheta;
    grid.nzeta = nzeta;
    grid.r.resize(static_cast<std::size_t>(ntheta * nzeta));
    grid.z.resize(grid.r.size());
    grid.nu.resize(grid.r.size());
    for (int zeta_index = 0; zeta_index < nzeta; ++zeta_index) {
        const double zeta = 2.0 * std::numbers::pi * zeta_index /
                            static_cast<double>(nzeta);
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index /
                                 static_cast<double>(ntheta);
            const std::size_t point = static_cast<std::size_t>(
                theta_index + ntheta * zeta_index);
            grid.r[point] = 3.0 + 0.4 * std::cos(2.0 * theta + zeta);
            grid.z[point] = 0.2 * std::sin(theta - 2.0 * zeta);
            grid.nu[point] = -0.1 * std::sin(zeta);
        }
    }
    const auto spectrum =
        magnetic_coordinate::analyze_mixed_grid_gpu(grid, 3, 3);
    const double two_pi_squared =
        2.0 * std::numbers::pi * std::numbers::pi;
    expect_near(spectrum.r[find_mode(spectrum, 0, 0)],
                {12.0 * std::numbers::pi * std::numbers::pi, 0.0}, 1.0e-11,
                "R zero-mode normalization");
    expect_near(spectrum.r[find_mode(spectrum, 2, 1)],
                {0.4 * two_pi_squared, 0.0}, 1.0e-11,
                "R signed mode");
    expect_near(spectrum.z[find_mode(spectrum, 1, -2)],
                {0.0, -0.2 * two_pi_squared}, 1.0e-11,
                "Z signed mode");
    expect_near(spectrum.nu[find_mode(spectrum, 0, 1)],
                {0.0, 0.1 * two_pi_squared}, 1.0e-11,
                "nu signed mode");
}

}  // namespace

int main() {
    int devices = 0;
    const cudaError_t status = cudaGetDeviceCount(&devices);
    if (status == cudaErrorNoDevice || status == cudaErrorInsufficientDriver ||
        devices == 0) {
        std::cout << "mixed_spectrum_test: SKIP (no CUDA device)\n";
        return 0;
    }
    if (status != cudaSuccess) {
        std::cerr << "mixed_spectrum_test: cudaGetDeviceCount: "
                  << cudaGetErrorString(status) << '\n';
        return 1;
    }
    try {
        test_integral_normalization();
    } catch (const std::exception& error) {
        std::cerr << "mixed_spectrum_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "mixed_spectrum_test: PASS\n";
    return 0;
}
