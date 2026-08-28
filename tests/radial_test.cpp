#include "magnetic_coordinate/radial.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expect_near(double actual,
                 double expected,
                 double tolerance,
                 const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": got " + std::to_string(actual) +
                                 ", expected " + std::to_string(expected));
    }
}

double cubic(double s, double point) {
    return 2.0 - 0.2 * s + 0.08 * s * s - 0.005 * s * s * s + point;
}

void test_two_point_interpolation() {
    constexpr int ns = 5;
    constexpr std::size_t points = 3;
    std::vector<double> half(static_cast<std::size_t>(ns - 1) * points);
    for (int j = 0; j < ns - 1; ++j) {
        const double s = static_cast<double>(j) + 0.5;
        for (std::size_t point = 0; point < points; ++point) {
            half[static_cast<std::size_t>(j) * points + point] =
                1.25 + 0.4 * s + static_cast<double>(point);
        }
    }

    const auto full = magnetic_coordinate::interpolate_half_to_full(
        half, ns, points,
        magnetic_coordinate::RadialInterpolationOrder::TWO_POINT);
    for (int j = 0; j < ns; ++j) {
        for (std::size_t point = 0; point < points; ++point) {
            expect_near(full[static_cast<std::size_t>(j) * points + point],
                        1.25 + 0.4 * static_cast<double>(j) +
                            static_cast<double>(point),
                        1.0e-14, "two-point linear interpolation");
        }
    }
}

void test_four_point_interpolation() {
    constexpr int ns = 8;
    constexpr std::size_t points = 2;
    std::vector<double> half(static_cast<std::size_t>(ns - 1) * points);
    for (int j = 0; j < ns - 1; ++j) {
        const double s = static_cast<double>(j) + 0.5;
        for (std::size_t point = 0; point < points; ++point) {
            half[static_cast<std::size_t>(j) * points + point] =
                cubic(s, static_cast<double>(point));
        }
    }

    const auto full = magnetic_coordinate::interpolate_half_to_full(
        half, ns, points,
        magnetic_coordinate::RadialInterpolationOrder::FOUR_POINT);
    for (int j = 0; j < ns; ++j) {
        for (std::size_t point = 0; point < points; ++point) {
            expect_near(
                full[static_cast<std::size_t>(j) * points + point],
                cubic(static_cast<double>(j), static_cast<double>(point)),
                2.0e-13, "four-point cubic interpolation");
        }
    }
}

void test_prepare_integer_grid() {
    constexpr int ns = 7;
    constexpr int ntheta = 4;
    constexpr int nzeta = 3;
    constexpr std::size_t points = static_cast<std::size_t>(ntheta * nzeta);
    const std::size_t half_size = static_cast<std::size_t>(ns - 1) * points;

    std::vector<double> sqrtg(half_size);
    std::vector<double> bsups(half_size, 0.0);
    std::vector<double> bsupu(half_size);
    std::vector<double> bsupv(half_size);
    std::vector<double> bsubs(half_size, 0.0);
    std::vector<double> bsubu(half_size, 0.0);
    std::vector<double> bsubv(half_size);

    for (int j = 0; j < ns - 1; ++j) {
        const double s = static_cast<double>(j) + 0.5;
        const double jacobian = -2.0 - 0.05 * s;
        const double phip = 1.5 + 0.1 * s;
        const double iota = 0.35 + 0.02 * s;
        const double chip = iota * phip;
        const double b2 = cubic(s, 0.0);
        for (std::size_t point = 0; point < points; ++point) {
            const std::size_t index =
                static_cast<std::size_t>(j) * points + point;
            sqrtg[index] = jacobian;
            bsupv[index] = phip / jacobian;
            bsupu[index] = chip / jacobian;
            bsubv[index] = b2 / bsupv[index];
        }
    }

    const magnetic_coordinate::NativeFieldView native{
        ns, ntheta, nzeta, sqrtg, bsups, bsupu, bsupv, bsubs, bsubu, bsubv,
    };
    const auto prepared = magnetic_coordinate::prepare_integer_grid(
        native, magnetic_coordinate::RadialInterpolationOrder::FOUR_POINT);

    check(prepared.b2.size() == static_cast<std::size_t>(ns) * points,
          "prepared B^2 extent");
    check(prepared.iota.size() == static_cast<std::size_t>(ns),
          "prepared iota extent");
    for (int j = 0; j < ns; ++j) {
        const double s = static_cast<double>(j);
        expect_near(
            prepared.toroidal_flux_derivative[static_cast<std::size_t>(j)],
            1.5 + 0.1 * s, 1.0e-13, "recovered toroidal flux derivative");
        expect_near(prepared.iota[static_cast<std::size_t>(j)], 0.35 + 0.02 * s,
                    2.0e-13, "recovered iota");
        for (std::size_t point = 0; point < points; ++point) {
            const std::size_t index =
                static_cast<std::size_t>(j) * points + point;
            expect_near(prepared.b2[index], cubic(s, 0.0), 2.0e-13,
                        "prepared B^2");
            expect_near(prepared.sqrtg[index], -2.0 - 0.05 * s, 1.0e-13,
                        "prepared signed Jacobian");
        }
    }
}

void test_rejects_malformed_field() {
    std::vector<double> values(8, 1.0);
    const magnetic_coordinate::NativeFieldView native{
        3,      2,
        2,      values,
        values, values,
        values, values,
        values, std::span<const double>(values).first(7),
    };

    bool rejected = false;
    try {
        static_cast<void>(magnetic_coordinate::prepare_integer_grid(
            native, magnetic_coordinate::RadialInterpolationOrder::TWO_POINT));
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "malformed native field must be rejected");
}

}  // namespace

int main() {
    try {
        test_two_point_interpolation();
        test_four_point_interpolation();
        test_prepare_integer_grid();
        test_rejects_malformed_field();
    } catch (const std::exception& error) {
        std::cerr << "radial_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "radial_test: PASS\n";
    return 0;
}
