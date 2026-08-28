#include "magnetic_coordinate/mixed_grid.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect_near(double actual, double expected, double tolerance,
                 const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": got " + std::to_string(actual) +
                                 ", expected " + std::to_string(expected));
    }
}

void test_analytic_double_remap() {
    constexpr int ns = 3;
    constexpr int source_ntheta = 128;
    constexpr int output_ntheta = 96;
    constexpr int nzeta = 2;
    constexpr double iota_value = 0.4;
    constexpr double lambda_amplitude = 0.08;
    constexpr double nu_amplitude = 0.12;
    constexpr double q_zero = -2.5;
    constexpr double four_pi_squared =
        4.0 * std::numbers::pi * std::numbers::pi;
    const std::size_t source_points =
        static_cast<std::size_t>(source_ntheta * nzeta);
    const std::size_t native_count = static_cast<std::size_t>(ns) * source_points;
    const std::size_t transformed_count =
        static_cast<std::size_t>(ns - 1) * source_points;

    magnetic_coordinate::NativeAngularGeometry geometry;
    geometry.ns = ns;
    geometry.ntheta = source_ntheta;
    geometry.nzeta = nzeta;
    geometry.r.resize(native_count);
    geometry.z.resize(native_count);
    geometry.lambda.resize(native_count);

    magnetic_coordinate::PestGrid pest;
    pest.source_ns = ns;
    pest.ntheta = source_ntheta;
    pest.nzeta = nzeta;
    pest.b2.resize(transformed_count);
    pest.sqrtg.resize(transformed_count);

    magnetic_coordinate::BoozerShift shift;
    shift.source_ns = ns;
    shift.ntheta = source_ntheta;
    shift.nzeta = nzeta;
    shift.nu.resize(transformed_count);
    shift.b2j00.assign(static_cast<std::size_t>(ns - 1),
                       four_pi_squared * q_zero);

    for (int surface = 0; surface < ns; ++surface) {
        for (int zeta = 0; zeta < nzeta; ++zeta) {
            for (int index = 0; index < source_ntheta; ++index) {
                const double theta = 2.0 * std::numbers::pi * index /
                                     static_cast<double>(source_ntheta);
                const std::size_t point = static_cast<std::size_t>(index) +
                    static_cast<std::size_t>(source_ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(nzeta * surface));
                geometry.lambda[point] = lambda_amplitude * std::sin(theta);
                geometry.r[point] = 10.0 + std::cos(theta);
                geometry.z[point] = std::sin(theta);
            }
        }
    }
    for (int surface = 0; surface < ns - 1; ++surface) {
        for (int zeta = 0; zeta < nzeta; ++zeta) {
            for (int index = 0; index < source_ntheta; ++index) {
                const double theta_p = 2.0 * std::numbers::pi * index /
                                       static_cast<double>(source_ntheta);
                const std::size_t point = static_cast<std::size_t>(index) +
                    static_cast<std::size_t>(source_ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(nzeta * surface));
                shift.nu[point] = nu_amplitude * std::sin(theta_p);
                pest.b2[point] = 4.0 + 0.2 * std::cos(theta_p);
                pest.sqrtg[point] = q_zero / pest.b2[point];
            }
        }
    }

    const std::vector<double> iota{0.0, iota_value, iota_value};
    const auto mixed = magnetic_coordinate::remap_to_boozer_mixed_grid(
        geometry, pest, shift, iota, output_ntheta);
    for (std::size_t point = 0; point < mixed.r.size(); ++point) {
        const double theta_p = mixed.source_theta_p[point];
        const double theta = mixed.source_theta[point];
        const int theta_b_index = static_cast<int>(point % output_ntheta);
        const double theta_b = 2.0 * std::numbers::pi * theta_b_index /
                               static_cast<double>(output_ntheta);
        expect_near(std::remainder(theta_p + iota_value * nu_amplitude *
                                      std::sin(theta_p) - theta_b,
                                  2.0 * std::numbers::pi),
                    0.0, 3.0e-7, "Boozer inverse root");
        expect_near(std::remainder(theta + lambda_amplitude * std::sin(theta) -
                                      theta_p,
                                  2.0 * std::numbers::pi),
                    0.0, 3.0e-7, "native inverse root");
        const double b2 = 4.0 + 0.2 * std::cos(theta_p);
        expect_near(mixed.r[point], 10.0 + std::cos(theta), 2.0e-7,
                    "mixed R");
        expect_near(mixed.z[point], std::sin(theta), 2.0e-7, "mixed Z");
        expect_near(mixed.b[point], std::sqrt(b2), 2.0e-7, "mixed B");
        expect_near(mixed.sqrtg_b[point], q_zero / b2, 2.0e-7,
                    "recovered Boozer Jacobian");
    }
}

}  // namespace

int main() {
    try {
        test_analytic_double_remap();
    } catch (const std::exception& error) {
        std::cerr << "mixed_grid_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "mixed_grid_test: PASS\n";
    return 0;
}
