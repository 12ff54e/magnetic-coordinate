#include "magnetic_coordinate/radial.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace magnetic_coordinate {
namespace {

std::size_t checked_product(std::size_t left,
                            std::size_t right,
                            const char* description) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::overflow_error(std::string(description) +
                                  " element count overflows size_t");
    }
    return left * right;
}

void require_finite(std::span<const double> values, const char* name) {
    for (double value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(std::string(name) +
                                        " contains a non-finite value");
        }
    }
}

double compensated_mean(std::span<const double> first,
                        std::span<const double> second) {
    if (first.size() != second.size() || first.empty()) {
        throw std::invalid_argument(
            "angular mean requires equal nonempty fields");
    }

    double sum = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < first.size(); ++i) {
        const double term = first[i] * second[i];
        const double adjusted = term - correction;
        const double next = sum + adjusted;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    return sum / static_cast<double>(first.size());
}

std::array<double, 4> interpolation_weights(int full_surface,
                                            int half_surfaces,
                                            int width,
                                            int& start) {
    const int preferred = full_surface - width / 2;
    start = std::clamp(preferred, 0, half_surfaces - width);

    std::array<double, 4> weights{};
    const double target = static_cast<double>(full_surface);
    for (int a = 0; a < width; ++a) {
        const double xa = static_cast<double>(start + a) + 0.5;
        double weight = 1.0;
        for (int b = 0; b < width; ++b) {
            if (a == b) continue;
            const double xb = static_cast<double>(start + b) + 0.5;
            weight *= (target - xb) / (xa - xb);
        }
        weights[static_cast<std::size_t>(a)] = weight;
    }
    return weights;
}

}  // namespace

std::vector<double> interpolate_half_to_full(
    std::span<const double> half_values,
    int ns,
    std::size_t values_per_surface,
    RadialInterpolationOrder order) {
    if (ns < 3) {
        throw std::invalid_argument(
            "radial interpolation requires at least three integer surfaces");
    }
    if (values_per_surface == 0) {
        throw std::invalid_argument(
            "radial interpolation requires a nonzero surface extent");
    }

    const int half_surfaces = ns - 1;
    const int width = static_cast<int>(order);
    if (width != 2 && width != 4) {
        throw std::invalid_argument("unsupported radial interpolation order");
    }
    if (half_surfaces < width) {
        throw std::invalid_argument(
            "radial grid is too small for the requested interpolation order");
    }

    const std::size_t expected =
        checked_product(static_cast<std::size_t>(half_surfaces),
                        values_per_surface, "half-grid radial");
    if (half_values.size() != expected) {
        throw std::invalid_argument("half-grid field has the wrong extent");
    }
    require_finite(half_values, "half-grid field");

    const std::size_t output_size = checked_product(
        static_cast<std::size_t>(ns), values_per_surface, "full-grid radial");
    std::vector<double> full_values(output_size, 0.0);
    for (int full_surface = 0; full_surface < ns; ++full_surface) {
        int start = 0;
        const auto weights =
            interpolation_weights(full_surface, half_surfaces, width, start);
        const std::size_t output_offset =
            static_cast<std::size_t>(full_surface) * values_per_surface;
        for (std::size_t point = 0; point < values_per_surface; ++point) {
            double value = 0.0;
            for (int stencil = 0; stencil < width; ++stencil) {
                const std::size_t input_offset =
                    static_cast<std::size_t>(start + stencil) *
                    values_per_surface;
                value += weights[static_cast<std::size_t>(stencil)] *
                         half_values[input_offset + point];
            }
            full_values[output_offset + point] = value;
        }
    }
    return full_values;
}

IntegerGridEquilibrium prepare_integer_grid(const NativeFieldView& native,
                                            RadialInterpolationOrder order) {
    if (native.ns < 3 || native.ntheta < 1 || native.nzeta < 1) {
        throw std::invalid_argument("invalid native equilibrium dimensions");
    }

    const std::size_t points =
        checked_product(static_cast<std::size_t>(native.ntheta),
                        static_cast<std::size_t>(native.nzeta), "angular");
    const std::size_t half_size = checked_product(
        static_cast<std::size_t>(native.ns - 1), points, "half-grid");

    const std::array<std::pair<std::span<const double>, const char*>, 7> fields{
        {
            {native.sqrtg, "sqrtg"},
            {native.bsups, "bsups"},
            {native.bsupu, "bsupu"},
            {native.bsupv, "bsupv"},
            {native.bsubs, "bsubs"},
            {native.bsubu, "bsubu"},
            {native.bsubv, "bsubv"},
        }};
    for (const auto& [field, name] : fields) {
        if (field.size() != half_size) {
            throw std::invalid_argument(std::string(name) +
                                        " has the wrong half-grid extent");
        }
        require_finite(field, name);
    }

    std::vector<double> b2_half(half_size, 0.0);
    for (std::size_t i = 0; i < half_size; ++i) {
        const double value = native.bsups[i] * native.bsubs[i] +
                             native.bsupu[i] * native.bsubu[i] +
                             native.bsupv[i] * native.bsubv[i];
        if (!(value > 0.0) || !std::isfinite(value)) {
            throw std::domain_error(
                "native field produced a nonpositive magnetic-field strength");
        }
        b2_half[i] = value;
    }

    std::vector<double> phip_half(static_cast<std::size_t>(native.ns - 1));
    std::vector<double> chip_half(static_cast<std::size_t>(native.ns - 1));
    for (int half_surface = 0; half_surface < native.ns - 1; ++half_surface) {
        const std::size_t offset =
            static_cast<std::size_t>(half_surface) * points;
        const auto sqrtg_surface = native.sqrtg.subspan(offset, points);
        phip_half[static_cast<std::size_t>(half_surface)] = compensated_mean(
            sqrtg_surface, native.bsupv.subspan(offset, points));
        chip_half[static_cast<std::size_t>(half_surface)] = compensated_mean(
            sqrtg_surface, native.bsupu.subspan(offset, points));
    }

    IntegerGridEquilibrium result;
    result.ns = native.ns;
    result.ntheta = native.ntheta;
    result.nzeta = native.nzeta;
    result.b2 = interpolate_half_to_full(b2_half, native.ns, points, order);
    result.sqrtg =
        interpolate_half_to_full(native.sqrtg, native.ns, points, order);
    result.toroidal_flux_derivative =
        interpolate_half_to_full(phip_half, native.ns, 1, order);
    result.poloidal_flux_derivative =
        interpolate_half_to_full(chip_half, native.ns, 1, order);

    result.b.resize(result.b2.size());
    for (std::size_t i = 0; i < result.b2.size(); ++i) {
        if (!(result.b2[i] > 0.0) || !std::isfinite(result.b2[i])) {
            throw std::domain_error(
                "radial interpolation produced a nonpositive magnetic-field "
                "strength");
        }
        result.b[i] = std::sqrt(result.b2[i]);
    }

    double flux_scale = 0.0;
    for (int surface = 0; surface < native.ns; ++surface) {
        const std::size_t index = static_cast<std::size_t>(surface);
        flux_scale = std::max(
            flux_scale,
            std::max(std::abs(result.toroidal_flux_derivative[index]),
                     std::abs(result.poloidal_flux_derivative[index])));
    }
    if (!(flux_scale > 0.0) || !std::isfinite(flux_scale)) {
        throw std::domain_error(
            "cannot recover iota from identically zero flux derivatives");
    }

    result.iota.resize(static_cast<std::size_t>(native.ns));
    const double flux_floor =
        64.0 * std::numeric_limits<double>::epsilon() * flux_scale;
    for (int surface = 0; surface < native.ns; ++surface) {
        const std::size_t index = static_cast<std::size_t>(surface);
        const double phip = result.toroidal_flux_derivative[index];
        const double chip = result.poloidal_flux_derivative[index];
        if (!std::isfinite(phip) || std::abs(phip) <= flux_floor) {
            throw std::domain_error(
                "cannot recover iota from a zero toroidal-flux derivative");
        }
        result.iota[index] = chip / phip;
    }
    return result;
}

}  // namespace magnetic_coordinate
