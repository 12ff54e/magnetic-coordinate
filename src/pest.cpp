#include "magnetic_coordinate/pest.hpp"

#include "magnetic_coordinate/detail/find_root.hpp"

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
    const double period = 2.0 * std::numbers::pi;
    double wrapped = std::fmod(angle, period);
    if (wrapped < 0.0) wrapped += period;
    return wrapped;
}

void validate_inputs(const NativeAngularGeometry& geometry,
                     const IntegerGridEquilibrium& fields) {
    if (geometry.ns != fields.ns || geometry.ntheta != fields.ntheta ||
        geometry.nzeta != fields.nzeta || geometry.ns < 2 ||
        geometry.ntheta < 4 || geometry.nzeta < 1) {
        throw std::invalid_argument("PEST input dimensions disagree");
    }
    const std::size_t points =
        static_cast<std::size_t>(geometry.ntheta) * geometry.nzeta;
    const std::size_t count = static_cast<std::size_t>(geometry.ns) * points;
    if (geometry.lambda.size() != count ||
        geometry.lambda_theta.size() != count || fields.b2.size() != count ||
        fields.sqrtg.size() != count) {
        throw std::invalid_argument("PEST input field extent mismatch");
    }
}

}  // namespace

PestGrid remap_to_pest(const NativeAngularGeometry& geometry,
                       const IntegerGridEquilibrium& integer_fields) {
    validate_inputs(geometry, integer_fields);
    const double period = 2.0 * std::numbers::pi;
    const std::size_t source_points =
        static_cast<std::size_t>(geometry.ntheta) * geometry.nzeta;
    const int transformed_surfaces = geometry.ns - 1;
    const std::size_t output_count =
        static_cast<std::size_t>(transformed_surfaces) * source_points;

    PestGrid result;
    result.source_ns = geometry.ns;
    result.ntheta = geometry.ntheta;
    result.nzeta = geometry.nzeta;
    result.source_theta.resize(output_count);
    result.b2.resize(output_count);
    result.sqrtg.resize(output_count);

    intp::InterpolationFunctionTemplate1D<SPLINE_ORDER> spline_template(
        std::make_pair(0.0, period), static_cast<std::size_t>(geometry.ntheta),
        true);
    std::vector<double> lambda_line(static_cast<std::size_t>(geometry.ntheta));
    std::vector<double> b2_line(static_cast<std::size_t>(geometry.ntheta));
    std::vector<double> sqrtg_line(static_cast<std::size_t>(geometry.ntheta));

    for (int surface = result.first_surface; surface < geometry.ns; ++surface) {
        const int output_surface = surface - result.first_surface;
        for (int zeta = 0; zeta < geometry.nzeta; ++zeta) {
            const std::size_t source_line =
                static_cast<std::size_t>(surface) * source_points +
                static_cast<std::size_t>(zeta) * geometry.ntheta;
            for (int theta = 0; theta < geometry.ntheta; ++theta) {
                const std::size_t line_index = static_cast<std::size_t>(theta);
                const std::size_t source = source_line + line_index;
                lambda_line[line_index] = geometry.lambda[source];
                b2_line[line_index] = integer_fields.b2[source];
                sqrtg_line[line_index] = integer_fields.sqrtg[source];
            }

            const auto lambda_spline = spline_template.interpolate(
                std::make_pair(lambda_line.begin(), lambda_line.end()));
            const auto b2_spline = spline_template.interpolate(
                std::make_pair(b2_line.begin(), b2_line.end()));
            const auto sqrtg_spline = spline_template.interpolate(
                std::make_pair(sqrtg_line.begin(), sqrtg_line.end()));

            const int validation_points = 4 * geometry.ntheta;
            for (int sample = 0; sample < validation_points; ++sample) {
                const double theta = period * static_cast<double>(sample) /
                                     static_cast<double>(validation_points);
                const double determinant =
                    1.0 + lambda_spline.derivative_at(
                              std::make_pair(theta, std::size_t{1}));
                if (!(determinant > 0.0) || !std::isfinite(determinant)) {
                    throw std::domain_error(
                        "theta-to-PEST map is not orientation preserving");
                }
            }

            for (int theta_p_index = 0; theta_p_index < geometry.ntheta;
                 ++theta_p_index) {
                const double theta_p = period *
                                       static_cast<double>(theta_p_index) /
                                       static_cast<double>(geometry.ntheta);
                const auto residual = [&](double theta) {
                    return theta + lambda_spline(theta) - theta_p;
                };
                const double root = detail::find_root(
                    residual, theta_p - period, theta_p + period,
                    [](double left, double right) {
                        return right - left <=
                               32.0 * std::numeric_limits<double>::epsilon() *
                                   (1.0 +
                                    std::max(std::abs(left), std::abs(right)));
                    });
                const double theta = wrap_angle(root);
                const double determinant =
                    1.0 + lambda_spline.derivative_at(
                              std::make_pair(theta, std::size_t{1}));
                if (!(determinant > 0.0) || !std::isfinite(determinant)) {
                    throw std::domain_error(
                        "invalid PEST Jacobian map determinant");
                }

                const std::size_t output =
                    static_cast<std::size_t>(theta_p_index) +
                    static_cast<std::size_t>(geometry.ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(geometry.nzeta) *
                             static_cast<std::size_t>(output_surface));
                result.source_theta[output] = theta;
                result.b2[output] = b2_spline(theta);
                result.sqrtg[output] = sqrtg_spline(theta) / determinant;
                if (!(result.b2[output] > 0.0) ||
                    !std::isfinite(result.b2[output]) ||
                    !std::isfinite(result.sqrtg[output])) {
                    throw std::domain_error(
                        "PEST interpolation produced a nonphysical field");
                }
            }
        }
    }
    return result;
}

}  // namespace magnetic_coordinate
