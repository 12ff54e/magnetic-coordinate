#include "magnetic_coordinate/angular_resample.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

#include <bsintp/Interpolation.hpp>

namespace magnetic_coordinate {
namespace {

constexpr std::size_t SPLINE_ORDER = 3;

void validate(const NativeAngularGeometry& geometry,
              const IntegerGridEquilibrium& fields,
              int output_nzeta) {
    if (geometry.ns != fields.ns || geometry.ntheta != fields.ntheta ||
        geometry.nzeta != fields.nzeta || geometry.ns < 2 ||
        geometry.ntheta < 1 || geometry.nzeta < 1 || output_nzeta < 1) {
        throw std::invalid_argument("toroidal resampling dimensions disagree");
    }
    if (geometry.nzeta != 1 && geometry.nzeta < 4) {
        throw std::invalid_argument(
            "periodic cubic toroidal resampling requires nzeta >= 4");
    }
    const std::size_t count = static_cast<std::size_t>(geometry.ns) *
                              geometry.ntheta * geometry.nzeta;
    const auto check = [count](const std::vector<double>& values) {
        if (values.size() != count) {
            throw std::invalid_argument(
                "toroidal resampling field extent mismatch");
        }
        for (double value : values) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "toroidal resampling field is non-finite");
            }
        }
    };
    check(geometry.r);
    check(geometry.z);
    check(geometry.lambda);
    check(geometry.lambda_theta);
    check(geometry.lambda_zeta);
    check(fields.b2);
    check(fields.b);
    check(fields.sqrtg);
}

std::vector<double> resample_field(const std::vector<double>& input,
                                   int ns,
                                   int ntheta,
                                   int source_nzeta,
                                   int output_nzeta) {
    if (source_nzeta == output_nzeta) return input;
    const std::size_t output_count =
        static_cast<std::size_t>(ns) * ntheta * output_nzeta;
    std::vector<double> output(output_count);
    if (source_nzeta == 1) {
        for (int surface = 0; surface < ns; ++surface) {
            for (int zeta = 0; zeta < output_nzeta; ++zeta) {
                for (int theta = 0; theta < ntheta; ++theta) {
                    const std::size_t source =
                        static_cast<std::size_t>(theta) +
                        static_cast<std::size_t>(ntheta) * surface;
                    const std::size_t target =
                        static_cast<std::size_t>(theta) +
                        static_cast<std::size_t>(ntheta) *
                            (static_cast<std::size_t>(zeta) +
                             static_cast<std::size_t>(output_nzeta) * surface);
                    output[target] = input[source];
                }
            }
        }
        return output;
    }

    constexpr double period = 2.0 * std::numbers::pi;
    intp::InterpolationFunctionTemplate1D<SPLINE_ORDER> spline_template(
        std::make_pair(0.0, period),
        static_cast<std::size_t>(source_nzeta), true);
    std::vector<double> line(static_cast<std::size_t>(source_nzeta));
    for (int surface = 0; surface < ns; ++surface) {
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int zeta = 0; zeta < source_nzeta; ++zeta) {
                const std::size_t source =
                    static_cast<std::size_t>(theta) +
                    static_cast<std::size_t>(ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(source_nzeta) * surface);
                line[static_cast<std::size_t>(zeta)] = input[source];
            }
            const auto spline = spline_template.interpolate(
                std::make_pair(line.begin(), line.end()));
            for (int zeta = 0; zeta < output_nzeta; ++zeta) {
                const double angle = period * static_cast<double>(zeta) /
                                     static_cast<double>(output_nzeta);
                const std::size_t target =
                    static_cast<std::size_t>(theta) +
                    static_cast<std::size_t>(ntheta) *
                        (static_cast<std::size_t>(zeta) +
                         static_cast<std::size_t>(output_nzeta) * surface);
                output[target] = spline(angle);
            }
        }
    }
    return output;
}

}  // namespace

ResampledAngularData resample_toroidal_grid(
    const NativeAngularGeometry& geometry,
    const IntegerGridEquilibrium& fields,
    int output_nzeta) {
    validate(geometry, fields, output_nzeta);
    ResampledAngularData result;
    result.geometry.ns = geometry.ns;
    result.geometry.ntheta = geometry.ntheta;
    result.geometry.nzeta = output_nzeta;
    result.geometry.r = resample_field(geometry.r, geometry.ns, geometry.ntheta,
                                       geometry.nzeta, output_nzeta);
    result.geometry.z = resample_field(geometry.z, geometry.ns, geometry.ntheta,
                                       geometry.nzeta, output_nzeta);
    result.geometry.lambda = resample_field(
        geometry.lambda, geometry.ns, geometry.ntheta, geometry.nzeta,
        output_nzeta);
    result.geometry.lambda_theta = resample_field(
        geometry.lambda_theta, geometry.ns, geometry.ntheta, geometry.nzeta,
        output_nzeta);
    result.geometry.lambda_zeta = resample_field(
        geometry.lambda_zeta, geometry.ns, geometry.ntheta, geometry.nzeta,
        output_nzeta);

    result.fields.ns = fields.ns;
    result.fields.ntheta = fields.ntheta;
    result.fields.nzeta = output_nzeta;
    result.fields.b2 = resample_field(fields.b2, fields.ns, fields.ntheta,
                                      fields.nzeta, output_nzeta);
    result.fields.sqrtg = resample_field(
        fields.sqrtg, fields.ns, fields.ntheta, fields.nzeta, output_nzeta);
    result.fields.b.resize(result.fields.b2.size());
    for (std::size_t index = 0; index < result.fields.b2.size(); ++index) {
        if (!(result.fields.b2[index] > 0.0)) {
            throw std::domain_error(
                "toroidal interpolation produced nonpositive B^2");
        }
        result.fields.b[index] = std::sqrt(result.fields.b2[index]);
    }
    result.fields.toroidal_flux_derivative =
        fields.toroidal_flux_derivative;
    result.fields.poloidal_flux_derivative =
        fields.poloidal_flux_derivative;
    result.fields.iota = fields.iota;
    return result;
}

}  // namespace magnetic_coordinate
