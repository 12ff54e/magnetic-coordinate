#include "magnetic_coordinate/pest.hpp"

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

void test_pest_remap() {
    constexpr int ns = 4;
    constexpr int ntheta = 128;
    constexpr int nzeta = 2;
    constexpr double amplitude = 0.2;
    constexpr double pest_jacobian = -3.0;
    const std::size_t points = static_cast<std::size_t>(ntheta * nzeta);
    const std::size_t count = static_cast<std::size_t>(ns) * points;

    magnetic_coordinate::NativeAngularGeometry geometry;
    geometry.ns = ns;
    geometry.ntheta = ntheta;
    geometry.nzeta = nzeta;
    geometry.lambda.resize(count);
    geometry.lambda_theta.resize(count);
    geometry.lambda_zeta.assign(count, 0.0);
    geometry.r.assign(count, 0.0);
    geometry.z.assign(count, 0.0);

    magnetic_coordinate::IntegerGridEquilibrium fields;
    fields.ns = ns;
    fields.ntheta = ntheta;
    fields.nzeta = nzeta;
    fields.b2.resize(count);
    fields.sqrtg.resize(count);

    for (int surface = 0; surface < ns; ++surface) {
        for (int zeta = 0; zeta < nzeta; ++zeta) {
            for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
                const double theta = 2.0 * std::numbers::pi * theta_index /
                                     static_cast<double>(ntheta);
                const std::size_t point =
                    static_cast<std::size_t>(theta_index) +
                    static_cast<std::size_t>(ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(nzeta * surface));
                geometry.lambda[point] = amplitude * std::sin(theta);
                geometry.lambda_theta[point] = amplitude * std::cos(theta);
                fields.b2[point] = 4.0 + 0.5 * std::sin(theta);
                fields.sqrtg[point] =
                    pest_jacobian * (1.0 + amplitude * std::cos(theta));
            }
        }
    }

    const auto pest = magnetic_coordinate::remap_to_pest(geometry, fields);
    for (int output_surface = 0; output_surface < ns - 1; ++output_surface) {
        for (int zeta = 0; zeta < nzeta; ++zeta) {
            for (int theta_p_index = 0; theta_p_index < ntheta;
                 ++theta_p_index) {
                const std::size_t point =
                    static_cast<std::size_t>(theta_p_index) +
                    static_cast<std::size_t>(ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(nzeta * output_surface));
                const double theta = pest.source_theta[point];
                const double theta_p = 2.0 * std::numbers::pi * theta_p_index /
                                       static_cast<double>(ntheta);
                double residual = theta + amplitude * std::sin(theta) - theta_p;
                residual = std::remainder(residual, 2.0 * std::numbers::pi);
                expect_near(residual, 0.0, 2.0e-7, "PEST inverse root");
                expect_near(pest.b2[point], 4.0 + 0.5 * std::sin(theta), 2.0e-7,
                            "PEST B^2 interpolation");
                expect_near(pest.sqrtg[point], pest_jacobian, 2.0e-6,
                            "PEST Jacobian transform");
            }
        }
    }
}

void test_rejects_folded_map() {
    constexpr int ns = 3;
    constexpr int ntheta = 32;
    constexpr int nzeta = 1;
    const std::size_t count = static_cast<std::size_t>(ns * ntheta * nzeta);
    magnetic_coordinate::NativeAngularGeometry geometry;
    geometry.ns = ns;
    geometry.ntheta = ntheta;
    geometry.nzeta = nzeta;
    geometry.lambda.resize(count);
    geometry.lambda_theta.resize(count);
    geometry.lambda_zeta.assign(count, 0.0);
    magnetic_coordinate::IntegerGridEquilibrium fields;
    fields.ns = ns;
    fields.ntheta = ntheta;
    fields.nzeta = nzeta;
    fields.b2.assign(count, 1.0);
    fields.sqrtg.assign(count, -1.0);
    for (int surface = 0; surface < ns; ++surface) {
        for (int theta_index = 0; theta_index < ntheta; ++theta_index) {
            const double theta = 2.0 * std::numbers::pi * theta_index /
                                 static_cast<double>(ntheta);
            const std::size_t point =
                static_cast<std::size_t>(surface * ntheta + theta_index);
            geometry.lambda[point] = 1.2 * std::sin(theta);
            geometry.lambda_theta[point] = 1.2 * std::cos(theta);
        }
    }

    bool rejected = false;
    try {
        static_cast<void>(magnetic_coordinate::remap_to_pest(geometry, fields));
    } catch (const std::domain_error&) { rejected = true; }
    if (!rejected) throw std::runtime_error("folded PEST map was accepted");
}

}  // namespace

int main() {
    try {
        test_pest_remap();
        test_rejects_folded_map();
    } catch (const std::exception& error) {
        std::cerr << "pest_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "pest_test: PASS\n";
    return 0;
}
