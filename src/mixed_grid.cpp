#include "magnetic_coordinate/mixed_grid.hpp"

#include "magnetic_coordinate/detail/find_root.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include <bsintp/Interpolation.hpp>

namespace magnetic_coordinate {
namespace {

constexpr std::size_t SPLINE_ORDER = 3;

double wrap_angle(double angle) {
    constexpr double period = 2.0 * std::numbers::pi;
    double wrapped = std::fmod(angle, period);
    if (wrapped < 0.0) wrapped += period;
    return wrapped;
}

template <typename Function>
double invert_periodic_map(const Function& displacement,
                           double target,
                           const char* description) {
    constexpr double period = 2.0 * std::numbers::pi;
    const auto residual = [&](double angle) {
        return angle + displacement(angle) - target;
    };
    try {
        return wrap_angle(detail::find_root(
            residual, target - period, target + period,
            [](double left, double right) {
                return right - left <=
                       32.0 * std::numeric_limits<double>::epsilon() *
                           (1.0 + std::max(std::abs(left), std::abs(right)));
            }));
    } catch (const std::exception&) {
        throw std::domain_error(description);
    }
}

void validate_inputs(const NativeAngularGeometry& geometry,
                     const PestGrid& pest,
                     const BoozerShift& shift,
                     std::span<const double> iota,
                     int output_ntheta) {
    if (geometry.ns != pest.source_ns || pest.source_ns != shift.source_ns ||
        pest.first_surface != shift.first_surface ||
        geometry.ntheta != pest.ntheta || pest.ntheta != shift.ntheta ||
        geometry.nzeta != pest.nzeta || pest.nzeta != shift.nzeta ||
        output_ntheta < 4 ||
        iota.size() != static_cast<std::size_t>(geometry.ns)) {
        throw std::invalid_argument("mixed-grid input dimensions disagree");
    }
    const std::size_t source_points =
        static_cast<std::size_t>(geometry.ntheta) * geometry.nzeta;
    const std::size_t native_count =
        static_cast<std::size_t>(geometry.ns) * source_points;
    const std::size_t transformed_count =
        static_cast<std::size_t>(geometry.ns - pest.first_surface) *
        source_points;
    if (geometry.r.size() != native_count || geometry.z.size() != native_count ||
        geometry.lambda.size() != native_count ||
        pest.b2.size() != transformed_count ||
        shift.nu.size() != transformed_count ||
        shift.b2j00.size() !=
            static_cast<std::size_t>(geometry.ns - pest.first_surface)) {
        throw std::invalid_argument("mixed-grid input field extent mismatch");
    }
}

}  // namespace

BoozerMixedGrid remap_to_boozer_mixed_grid(
    const NativeAngularGeometry& geometry,
    const PestGrid& pest,
    const BoozerShift& shift,
    std::span<const double> iota,
    int output_ntheta) {
    validate_inputs(geometry, pest, shift, iota, output_ntheta);
    constexpr double period = 2.0 * std::numbers::pi;
    constexpr double four_pi_squared =
        4.0 * std::numbers::pi * std::numbers::pi;
    const int surfaces = geometry.ns - pest.first_surface;
    const std::size_t source_points =
        static_cast<std::size_t>(geometry.ntheta) * geometry.nzeta;
    const std::size_t output_points =
        static_cast<std::size_t>(output_ntheta) * geometry.nzeta;
    const std::size_t output_count =
        static_cast<std::size_t>(surfaces) * output_points;

    BoozerMixedGrid result;
    result.source_ns = geometry.ns;
    result.first_surface = pest.first_surface;
    result.ntheta = output_ntheta;
    result.nzeta = geometry.nzeta;
    result.source_theta_p.resize(output_count);
    result.source_theta.resize(output_count);
    result.r.resize(output_count);
    result.z.resize(output_count);
    result.nu.resize(output_count);
    result.b.resize(output_count);
    result.sqrtg_b.resize(output_count);
    result.b2j00 = shift.b2j00;

    intp::InterpolationFunctionTemplate1D<SPLINE_ORDER> native_template(
        std::make_pair(0.0, period),
        static_cast<std::size_t>(geometry.ntheta), true);
    std::vector<double> lambda_line(static_cast<std::size_t>(geometry.ntheta));
    std::vector<double> r_line(lambda_line.size());
    std::vector<double> z_line(lambda_line.size());
    std::vector<double> nu_line(lambda_line.size());
    std::vector<double> b2_line(lambda_line.size());

    for (int output_surface = 0; output_surface < surfaces; ++output_surface) {
        const int source_surface = output_surface + pest.first_surface;
        const double surface_iota = iota[static_cast<std::size_t>(source_surface)];
        if (!std::isfinite(surface_iota)) {
            throw std::invalid_argument("mixed-grid iota is non-finite");
        }
        const double zero_mode =
            shift.b2j00[static_cast<std::size_t>(output_surface)];
        if (!std::isfinite(zero_mode) || zero_mode == 0.0) {
            throw std::domain_error("mixed-grid b2j00 is zero or non-finite");
        }

        for (int zeta = 0; zeta < geometry.nzeta; ++zeta) {
            const std::size_t native_line =
                static_cast<std::size_t>(source_surface) * source_points +
                static_cast<std::size_t>(zeta) * geometry.ntheta;
            const std::size_t pest_line =
                static_cast<std::size_t>(output_surface) * source_points +
                static_cast<std::size_t>(zeta) * geometry.ntheta;
            for (int theta = 0; theta < geometry.ntheta; ++theta) {
                const std::size_t local = static_cast<std::size_t>(theta);
                lambda_line[local] = geometry.lambda[native_line + local];
                r_line[local] = geometry.r[native_line + local];
                z_line[local] = geometry.z[native_line + local];
                nu_line[local] = shift.nu[pest_line + local];
                b2_line[local] = pest.b2[pest_line + local];
            }

            const auto lambda_spline = native_template.interpolate(
                std::make_pair(lambda_line.begin(), lambda_line.end()));
            const auto r_spline = native_template.interpolate(
                std::make_pair(r_line.begin(), r_line.end()));
            const auto z_spline = native_template.interpolate(
                std::make_pair(z_line.begin(), z_line.end()));
            const auto nu_spline = native_template.interpolate(
                std::make_pair(nu_line.begin(), nu_line.end()));
            const auto b2_spline = native_template.interpolate(
                std::make_pair(b2_line.begin(), b2_line.end()));

            const int validation_points = 4 * std::max(geometry.ntheta,
                                                        output_ntheta);
            for (int sample = 0; sample < validation_points; ++sample) {
                const double theta_p = period * static_cast<double>(sample) /
                                       static_cast<double>(validation_points);
                const double determinant =
                    1.0 + surface_iota * nu_spline.derivative_at(
                        std::make_pair(theta_p, std::size_t{1}));
                if (!(determinant > 0.0) || !std::isfinite(determinant)) {
                    throw std::domain_error(
                        "PEST-to-Boozer map is not orientation preserving");
                }
            }

            for (int theta_b_index = 0; theta_b_index < output_ntheta;
                 ++theta_b_index) {
                const double theta_b =
                    period * static_cast<double>(theta_b_index) /
                    static_cast<double>(output_ntheta);
                const auto iota_nu = [&](double theta_p) {
                    return surface_iota * nu_spline(theta_p);
                };
                const double theta_p = invert_periodic_map(
                    iota_nu, theta_b, "could not invert PEST-to-Boozer map");
                const double theta = invert_periodic_map(
                    lambda_spline, theta_p,
                    "could not recover the original poloidal angle");
                const double b2 = b2_spline(theta_p);
                if (!(b2 > 0.0) || !std::isfinite(b2)) {
                    throw std::domain_error(
                        "mixed-grid interpolation produced invalid B^2");
                }
                const std::size_t output =
                    static_cast<std::size_t>(theta_b_index) +
                    static_cast<std::size_t>(output_ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(geometry.nzeta) *
                             static_cast<std::size_t>(output_surface));
                result.source_theta_p[output] = theta_p;
                result.source_theta[output] = theta;
                result.r[output] = r_spline(theta);
                result.z[output] = z_spline(theta);
                result.nu[output] = nu_spline(theta_p);
                result.b[output] = std::sqrt(b2);
                result.sqrtg_b[output] = zero_mode / (four_pi_squared * b2);
                if (!std::isfinite(result.r[output]) ||
                    !std::isfinite(result.z[output]) ||
                    !std::isfinite(result.nu[output]) ||
                    !std::isfinite(result.sqrtg_b[output])) {
                    throw std::domain_error(
                        "mixed-grid interpolation produced non-finite output");
                }
            }
        }
    }
    return result;
}

}  // namespace magnetic_coordinate
