#include "magnetic_coordinate/spectral.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>

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

magnetic_coordinate::CumesEquilibrium make_equilibrium() {
    magnetic_coordinate::CumesEquilibrium equilibrium;
    equilibrium.ns = 5;
    equilibrium.mpol = 2;
    equilibrium.ntor = 1;
    equilibrium.nfp = 5;
    equilibrium.mnmax = equilibrium.mpol * (equilibrium.ntor + 1);
    equilibrium.ntheta = 8;
    equilibrium.nzeta = 6;
    equilibrium.phiedge = -2.5;
    equilibrium.aphi = {1.0, 0.25};
    for (auto& family : equilibrium.families) {
        family.assign(
            static_cast<std::size_t>(equilibrium.ns * equilibrium.mnmax), 0.0);
    }
    return equilibrium;
}

void set_coefficient(
    magnetic_coordinate::CumesEquilibrium& equilibrium,
    magnetic_coordinate::CumesEquilibrium::SpectralFamily family,
    int m,
    int n,
    int surface,
    double value) {
    const int mode = m * (equilibrium.ntor + 1) + n;
    equilibrium.families[family][static_cast<std::size_t>(
        surface + mode * equilibrium.ns)] = value;
}

void test_flux_normalization() {
    const auto equilibrium = make_equilibrium();
    const auto normalization =
        magnetic_coordinate::reconstruct_flux_normalization(equilibrium);

    const double maximum_flux = 1.0 / std::numbers::pi;
    expect_near(normalization.maximum_toroidal_flux, maximum_flux, 1.0e-15,
                "maximum toroidal flux");
    for (int surface = 0; surface < equilibrium.ns; ++surface) {
        const double s = static_cast<double>(surface) /
                         static_cast<double>(equilibrium.ns - 1);
        expect_near(normalization.phip_full[static_cast<std::size_t>(surface)],
                    maximum_flux * (1.0 + 0.5 * s), 1.0e-15,
                    "full-grid Phi-prime");
    }

    double expected_norm_squared = 0.0;
    for (int half_surface = 0; half_surface < equilibrium.ns - 1;
         ++half_surface) {
        const double s = (static_cast<double>(half_surface) + 0.5) /
                         static_cast<double>(equilibrium.ns - 1);
        const double phip = maximum_flux * (1.0 + 0.5 * s);
        expected_norm_squared +=
            phip * phip / static_cast<double>(equilibrium.ns - 1);
    }
    expect_near(normalization.lambda_scale, std::sqrt(expected_norm_squared),
                1.0e-15, "lambda scale");
}

void test_single_mode_synthesis() {
    auto equilibrium = make_equilibrium();
    constexpr int surface = 2;
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::RMNCC,
                    0, 0, surface, 10.0);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::RMNCC,
                    1, 1, surface, 0.5);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::RMNSS,
                    1, 1, surface, 0.2);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::ZMNSC,
                    1, 1, surface, -0.7);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::ZMNCS,
                    1, 1, surface, 0.4);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::LMNSC,
                    1, 1, surface, 0.1);
    set_coefficient(equilibrium, magnetic_coordinate::CumesEquilibrium::LMNCS,
                    1, 1, surface, -0.3);

    const auto normalization =
        magnetic_coordinate::reconstruct_flux_normalization(equilibrium);
    const auto geometry = magnetic_coordinate::synthesize_native_geometry(
        equilibrium, normalization);

    constexpr int theta_index = 1;
    constexpr int zeta_index = 2;
    const double theta = 2.0 * std::numbers::pi * theta_index /
                         static_cast<double>(equilibrium.ntheta);
    const double zeta = 2.0 * std::numbers::pi * zeta_index /
                        static_cast<double>(equilibrium.nzeta);
    const double cos_m = std::cos(theta);
    const double sin_m = std::sin(theta);
    const double cos_n = std::cos(zeta);
    const double sin_n = std::sin(zeta);
    const double factor =
        normalization.lambda_multiplier[static_cast<std::size_t>(surface)];
    const std::size_t point =
        static_cast<std::size_t>(theta_index) +
        static_cast<std::size_t>(equilibrium.ntheta) *
            (static_cast<std::size_t>(zeta_index) +
             static_cast<std::size_t>(equilibrium.nzeta * surface));

    expect_near(geometry.r[point],
                10.0 + 0.5 * cos_m * cos_n + 0.2 * sin_m * sin_n, 1.0e-14,
                "R synthesis");
    expect_near(geometry.z[point], -0.7 * sin_m * cos_n + 0.4 * cos_m * sin_n,
                1.0e-14, "Z synthesis");
    expect_near(geometry.lambda[point],
                factor * (0.1 * sin_m * cos_n - 0.3 * cos_m * sin_n), 1.0e-14,
                "lambda synthesis");
    expect_near(geometry.lambda_theta[point],
                factor * (0.1 * cos_m * cos_n + 0.3 * sin_m * sin_n), 1.0e-14,
                "lambda theta derivative");
    expect_near(geometry.lambda_zeta[point],
                factor * 5.0 * (-0.1 * sin_m * sin_n - 0.3 * cos_m * cos_n),
                1.0e-14, "lambda physical zeta derivative");
}

}  // namespace

int main() {
    try {
        test_flux_normalization();
        test_single_mode_synthesis();
    } catch (const std::exception& error) {
        std::cerr << "spectral_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "spectral_test: PASS\n";
    return 0;
}
