#include "magnetic_coordinate/boozer_fft.hpp"

#include <cuda_runtime.h>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect_near(double actual,
                 double expected,
                 double tolerance,
                 const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": got " + std::to_string(actual) +
                                 ", expected " + std::to_string(expected));
    }
}

void test_manufactured_shift() {
    constexpr int ntheta = 32;
    constexpr int nzeta = 16;
    constexpr int nfp = 5;
    constexpr double iota = 0.4;
    constexpr double constant = -3.0;
    constexpr double poloidal_amplitude = 0.12;
    constexpr double toroidal_amplitude = -0.08;
    magnetic_coordinate::PestGrid pest;
    pest.source_ns = 2;
    pest.first_surface = 1;
    pest.ntheta = ntheta;
    pest.nzeta = nzeta;
    pest.b2.assign(static_cast<std::size_t>(ntheta * nzeta), 1.0);
    pest.sqrtg.resize(pest.b2.size());
    for (int zeta_index = 0; zeta_index < nzeta; ++zeta_index) {
        const double zeta =
            2.0 * std::numbers::pi * zeta_index / static_cast<double>(nzeta);
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index /
                                 static_cast<double>(ntheta);
            const std::size_t point =
                static_cast<std::size_t>(theta_index + ntheta * zeta_index);
            pest.sqrtg[point] =
                constant * (1.0 + poloidal_amplitude * std::cos(2.0 * theta) +
                            toroidal_amplitude * std::cos(zeta));
        }
    }

    const std::vector<double> iota_profile{0.0, iota};
    const auto shift =
        magnetic_coordinate::solve_boozer_shift_gpu(pest, iota_profile, nfp);
    expect_near(shift.b2j00[0],
                4.0 * std::numbers::pi * std::numbers::pi * constant, 1.0e-11,
                "b2j00 normalization");
    for (int zeta_index = 0; zeta_index < nzeta; ++zeta_index) {
        const double zeta =
            2.0 * std::numbers::pi * zeta_index / static_cast<double>(nzeta);
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index /
                                 static_cast<double>(ntheta);
            const std::size_t point =
                static_cast<std::size_t>(theta_index + ntheta * zeta_index);
            const double expected =
                poloidal_amplitude / (2.0 * iota) * std::sin(2.0 * theta) +
                toroidal_amplitude / static_cast<double>(nfp) * std::sin(zeta);
            expect_near(shift.nu[point], expected, 2.0e-13,
                        "manufactured Boozer shift");
        }
    }
}

}  // namespace

int main() {
    int devices = 0;
    const cudaError_t status = cudaGetDeviceCount(&devices);
    if (status == cudaErrorNoDevice || status == cudaErrorInsufficientDriver ||
        devices == 0) {
        std::cout << "boozer_fft_test: SKIP (no CUDA device)\n";
        return 0;
    }
    if (status != cudaSuccess) {
        std::cerr << "boozer_fft_test: cudaGetDeviceCount: "
                  << cudaGetErrorString(status) << '\n';
        return 1;
    }
    try {
        test_manufactured_shift();
    } catch (const std::exception& error) {
        std::cerr << "boozer_fft_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "boozer_fft_test: PASS\n";
    return 0;
}
