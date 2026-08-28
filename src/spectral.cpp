#include "magnetic_coordinate/spectral.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>

namespace magnetic_coordinate {
namespace {

double toroidal_flux_derivative(std::span<const double> coefficients,
                                double s) {
    double derivative = 0.0;
    double power = 1.0;
    for (std::size_t degree = 0; degree < coefficients.size(); ++degree) {
        derivative +=
            static_cast<double>(degree + 1) * coefficients[degree] * power;
        power *= s;
    }
    return derivative;
}

double toroidal_flux_at_one(std::span<const double> coefficients) {
    double value = 0.0;
    for (double coefficient : coefficients) value += coefficient;
    return value;
}

void validate_spectral_input(const CumesEquilibriumView& equilibrium) {
    if (equilibrium.ns < 3 || equilibrium.mpol < 1 || equilibrium.ntor < 0 ||
        equilibrium.nfp < 1 || equilibrium.ntheta < 1 ||
        equilibrium.nzeta < 1 || equilibrium.mnmax < 1) {
        throw std::invalid_argument("invalid spectral equilibrium dimensions");
    }
    const std::size_t modes = static_cast<std::size_t>(equilibrium.mpol) *
                              static_cast<std::size_t>(equilibrium.ntor + 1);
    if (modes != static_cast<std::size_t>(equilibrium.mnmax)) {
        throw std::invalid_argument("spectral mode dimensions disagree");
    }
    const std::size_t family_size =
        modes * static_cast<std::size_t>(equilibrium.ns);
    for (const auto& family : equilibrium.families) {
        if (family.size() != family_size) {
            throw std::invalid_argument("spectral family has the wrong extent");
        }
        for (double value : family) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "spectral family contains a non-finite value");
            }
        }
    }
}

}  // namespace

FluxNormalization reconstruct_flux_normalization(
    const CumesEquilibriumView& equilibrium) {
    if (equilibrium.ns < 3 || equilibrium.aphi.empty() ||
        !std::isfinite(equilibrium.phiedge)) {
        throw std::invalid_argument(
            "cannot reconstruct flux normalization from incomplete metadata");
    }
    for (double coefficient : equilibrium.aphi) {
        if (!std::isfinite(coefficient)) {
            throw std::invalid_argument(
                "toroidal flux coefficients contain a non-finite value");
        }
    }

    const double flux_at_one = toroidal_flux_at_one(equilibrium.aphi);
    const auto largest_coefficient =
        std::max_element(equilibrium.aphi.begin(), equilibrium.aphi.end(),
                         [](double left, double right) {
                             return std::abs(left) < std::abs(right);
                         });
    const double profile_scale =
        std::max(std::abs(equilibrium.phiedge), std::abs(*largest_coefficient));
    const double floor =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(std::abs(profile_scale), std::numeric_limits<double>::min());
    if (std::abs(flux_at_one) <= floor) {
        throw std::domain_error("toroidal flux profile vanishes at the edge");
    }

    FluxNormalization result;
    result.maximum_toroidal_flux =
        -equilibrium.phiedge / (2.0 * std::numbers::pi * flux_at_one);
    const double delta_s = 1.0 / static_cast<double>(equilibrium.ns - 1);
    result.phip_full.resize(static_cast<std::size_t>(equilibrium.ns));
    result.phip_half.resize(static_cast<std::size_t>(equilibrium.ns - 1));

    for (int surface = 0; surface < equilibrium.ns; ++surface) {
        const double s = delta_s * static_cast<double>(surface);
        result.phip_full[static_cast<std::size_t>(surface)] =
            result.maximum_toroidal_flux *
            toroidal_flux_derivative(equilibrium.aphi, s);
    }
    double lambda_norm_squared = 0.0;
    for (int half_surface = 0; half_surface < equilibrium.ns - 1;
         ++half_surface) {
        const double s = delta_s * (static_cast<double>(half_surface) + 0.5);
        const double phip = result.maximum_toroidal_flux *
                            toroidal_flux_derivative(equilibrium.aphi, s);
        result.phip_half[static_cast<std::size_t>(half_surface)] = phip;
        lambda_norm_squared += phip * phip;
    }
    result.lambda_scale = std::sqrt(delta_s * lambda_norm_squared);
    if (!(result.lambda_scale > 0.0) || !std::isfinite(result.lambda_scale)) {
        throw std::domain_error("invalid cuMES lambda normalization");
    }

    double phip_scale = 0.0;
    for (double phip : result.phip_full) {
        phip_scale = std::max(phip_scale, std::abs(phip));
    }
    const double phip_floor =
        64.0 * std::numeric_limits<double>::epsilon() * phip_scale;
    result.lambda_multiplier.resize(static_cast<std::size_t>(equilibrium.ns));
    for (int surface = 0; surface < equilibrium.ns; ++surface) {
        const std::size_t index = static_cast<std::size_t>(surface);
        const double phip = result.phip_full[index];
        if (std::abs(phip) <= phip_floor || !std::isfinite(phip)) {
            throw std::domain_error(
                "physical PEST lambda is undefined where Phi-prime vanishes");
        }
        result.lambda_multiplier[index] = result.lambda_scale / phip;
    }
    return result;
}

NativeAngularGeometry synthesize_native_geometry(
    const CumesEquilibriumView& equilibrium,
    const FluxNormalization& normalization) {
    validate_spectral_input(equilibrium);
    if (normalization.lambda_multiplier.size() !=
        static_cast<std::size_t>(equilibrium.ns)) {
        throw std::invalid_argument("lambda normalization extent mismatch");
    }

    const std::size_t points = static_cast<std::size_t>(equilibrium.ntheta) *
                               static_cast<std::size_t>(equilibrium.nzeta);
    const std::size_t count = static_cast<std::size_t>(equilibrium.ns) * points;
    NativeAngularGeometry result;
    result.ns = equilibrium.ns;
    result.ntheta = equilibrium.ntheta;
    result.nzeta = equilibrium.nzeta;
    result.r.assign(count, 0.0);
    result.z.assign(count, 0.0);
    result.lambda.assign(count, 0.0);
    result.lambda_theta.assign(count, 0.0);
    result.lambda_zeta.assign(count, 0.0);

    const double two_pi = 2.0 * std::numbers::pi;
    for (int surface = 0; surface < equilibrium.ns; ++surface) {
        const double lambda_multiplier =
            normalization.lambda_multiplier[static_cast<std::size_t>(surface)];
        if (!std::isfinite(lambda_multiplier)) {
            throw std::invalid_argument(
                "lambda normalization contains a non-finite value");
        }
        for (int zeta_index = 0; zeta_index < equilibrium.nzeta; ++zeta_index) {
            const double zeta = two_pi * static_cast<double>(zeta_index) /
                                static_cast<double>(equilibrium.nzeta);
            for (int theta_index = 0; theta_index < equilibrium.ntheta;
                 ++theta_index) {
                const double theta = two_pi * static_cast<double>(theta_index) /
                                     static_cast<double>(equilibrium.ntheta);
                const std::size_t point =
                    static_cast<std::size_t>(theta_index) +
                    static_cast<std::size_t>(equilibrium.ntheta) *
                        (static_cast<std::size_t>(zeta_index) +
                         static_cast<std::size_t>(equilibrium.nzeta) *
                             static_cast<std::size_t>(surface));

                for (int m = 0; m < equilibrium.mpol; ++m) {
                    const double m_theta = static_cast<double>(m) * theta;
                    const double cos_m = std::cos(m_theta);
                    const double sin_m = std::sin(m_theta);
                    for (int n = 0; n <= equilibrium.ntor; ++n) {
                        const double n_zeta = static_cast<double>(n) * zeta;
                        const double cos_n = std::cos(n_zeta);
                        const double sin_n = std::sin(n_zeta);
                        const int mode = m * (equilibrium.ntor + 1) + n;
                        const std::size_t coefficient =
                            static_cast<std::size_t>(surface) +
                            static_cast<std::size_t>(mode) *
                                static_cast<std::size_t>(equilibrium.ns);

                        const double rcc =
                            equilibrium
                                .families[CumesEquilibrium::RMNCC][coefficient];
                        const double rss =
                            equilibrium
                                .families[CumesEquilibrium::RMNSS][coefficient];
                        const double zsc =
                            equilibrium
                                .families[CumesEquilibrium::ZMNSC][coefficient];
                        const double zcs =
                            equilibrium
                                .families[CumesEquilibrium::ZMNCS][coefficient];
                        const double lsc =
                            equilibrium
                                .families[CumesEquilibrium::LMNSC][coefficient];
                        const double lcs =
                            equilibrium
                                .families[CumesEquilibrium::LMNCS][coefficient];

                        result.r[point] +=
                            rcc * cos_m * cos_n + rss * sin_m * sin_n;
                        result.z[point] +=
                            zsc * sin_m * cos_n + zcs * cos_m * sin_n;
                        result.lambda[point] +=
                            lambda_multiplier *
                            (lsc * sin_m * cos_n + lcs * cos_m * sin_n);
                        result.lambda_theta[point] +=
                            lambda_multiplier * static_cast<double>(m) *
                            (lsc * cos_m * cos_n - lcs * sin_m * sin_n);
                        result.lambda_zeta[point] +=
                            lambda_multiplier * static_cast<double>(n) *
                            static_cast<double>(equilibrium.nfp) *
                            (-lsc * sin_m * sin_n + lcs * cos_m * cos_n);
                    }
                }
            }
        }
    }
    return result;
}

}  // namespace magnetic_coordinate
