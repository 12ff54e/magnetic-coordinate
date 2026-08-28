#include "magnetic_coordinate/angular_resample.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

void test_periodic_single_mode() {
    constexpr int ns = 3;
    constexpr int ntheta = 2;
    constexpr int source_nzeta = 16;
    constexpr int output_nzeta = 11;
    const std::size_t count =
        static_cast<std::size_t>(ns * ntheta * source_nzeta);
    magnetic_coordinate::NativeAngularGeometry geometry;
    geometry.ns = ns;
    geometry.ntheta = ntheta;
    geometry.nzeta = source_nzeta;
    geometry.r.resize(count);
    geometry.z.resize(count);
    geometry.lambda.resize(count);
    geometry.lambda_theta.resize(count);
    geometry.lambda_zeta.resize(count);
    magnetic_coordinate::IntegerGridEquilibrium fields;
    fields.ns = ns;
    fields.ntheta = ntheta;
    fields.nzeta = source_nzeta;
    fields.b2.resize(count);
    fields.b.resize(count);
    fields.sqrtg.resize(count);
    fields.iota.assign(ns, 0.4);
    fields.toroidal_flux_derivative.assign(ns, 1.0);
    fields.poloidal_flux_derivative.assign(ns, 0.4);
    for (int surface = 0; surface < ns; ++surface) {
        for (int zeta = 0; zeta < source_nzeta; ++zeta) {
            const double angle = 2.0 * std::numbers::pi * zeta /
                                 static_cast<double>(source_nzeta);
            for (int theta = 0; theta < ntheta; ++theta) {
                const std::size_t point = static_cast<std::size_t>(theta) +
                    static_cast<std::size_t>(ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(source_nzeta * surface));
                geometry.r[point] = 10.0 + 0.2 * std::cos(angle);
                geometry.z[point] = 0.3 * std::sin(angle);
                geometry.lambda[point] = 0.01 * std::sin(angle);
                geometry.lambda_theta[point] = 0.0;
                geometry.lambda_zeta[point] = 0.01 * std::cos(angle);
                fields.b2[point] = 4.0 + 0.1 * std::cos(angle);
                fields.b[point] = std::sqrt(fields.b2[point]);
                fields.sqrtg[point] = -2.0 + 0.1 * std::sin(angle);
            }
        }
    }
    const auto output = magnetic_coordinate::resample_toroidal_grid(
        geometry, fields, output_nzeta);
    for (int zeta = 0; zeta < output_nzeta; ++zeta) {
        const double angle = 2.0 * std::numbers::pi * zeta /
                             static_cast<double>(output_nzeta);
        const std::size_t point = static_cast<std::size_t>(ntheta * zeta);
        if (std::abs(output.geometry.r[point] -
                     (10.0 + 0.2 * std::cos(angle))) > 2.0e-5 ||
            std::abs(output.geometry.z[point] - 0.3 * std::sin(angle)) >
                2.0e-5 ||
            std::abs(output.fields.b2[point] -
                     (4.0 + 0.1 * std::cos(angle))) > 2.0e-5) {
            throw std::runtime_error("periodic toroidal resampling mismatch");
        }
    }
}

}  // namespace

int main() {
    try {
        test_periodic_single_mode();
    } catch (const std::exception& error) {
        std::cerr << "angular_resample_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "angular_resample_test: PASS\n";
    return 0;
}
