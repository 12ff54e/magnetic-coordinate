#include "magnetic_coordinate/boozer_output.hpp"

#include "boozer_output_internal.hpp"
#include "magnetic_coordinate/boozer_binary.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>

namespace magnetic_coordinate {
namespace {

std::size_t checked_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::invalid_argument("Boozer result dimension overflow");
    }
    return left * right;
}

std::string lower_extension(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char value) {
                       return static_cast<char>(std::tolower(value));
                   });
    return extension;
}

}  // namespace

BoozerOutputSpec resolve_boozer_output_spec(const std::filesystem::path& path) {
    const std::string extension = lower_extension(path);
    if (extension == ".bin") { return {BoozerOutputFormat::BINARY, path}; }
    if (extension == ".nc") { return {BoozerOutputFormat::NETCDF, path}; }
    if (extension == ".h5" || extension == ".hdf5") {
        return {BoozerOutputFormat::HDF5, path};
    }
    throw std::invalid_argument("unrecognized Boozer output suffix '" +
                                extension +
                                "' (expected .bin, .nc, .h5, or .hdf5)");
}

bool boozer_output_format_available(BoozerOutputFormat format) {
    switch (format) {
        case BoozerOutputFormat::BINARY:
            return true;
        case BoozerOutputFormat::NETCDF:
#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
            return true;
#else
            return false;
#endif
        case BoozerOutputFormat::HDF5:
#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
            return true;
#else
            return false;
#endif
    }
    return false;
}

std::string_view boozer_output_suffix(BoozerOutputFormat format) {
    switch (format) {
        case BoozerOutputFormat::BINARY:
            return ".bin";
        case BoozerOutputFormat::NETCDF:
            return ".nc";
        case BoozerOutputFormat::HDF5:
            return ".h5";
    }
    return "";
}

void write_boozer_output(const BoozerOutputSpec& spec,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    if (!boozer_output_format_available(spec.format)) {
        throw std::runtime_error(
            "Boozer output format '" +
            std::string(boozer_output_suffix(spec.format)) +
            "' is not available in this build");
    }
    switch (spec.format) {
        case BoozerOutputFormat::BINARY:
            write_boozer_binary(spec.path, result, source_path);
            return;
        case BoozerOutputFormat::NETCDF:
#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
            detail::write_boozer_netcdf(spec.path, result, source_path);
            return;
#else
            break;
#endif
        case BoozerOutputFormat::HDF5:
#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
            detail::write_boozer_hdf5(spec.path, result, source_path);
            return;
#else
            break;
#endif
    }
    throw std::runtime_error("unavailable Boozer output format");
}

namespace detail {

void validate_boozer_result(const BoozerResult& result) {
    const int surfaces = result.source_ns - result.grid.first_surface;
    if (result.spectrum.mmax < 0 || result.spectrum.nmax < 0 || surfaces < 1) {
        throw std::invalid_argument("inconsistent Boozer result dimensions");
    }
    const std::size_t mode_count =
        checked_product(2 * static_cast<std::size_t>(result.spectrum.mmax) + 1,
                        2 * static_cast<std::size_t>(result.spectrum.nmax) + 1);
    if (mode_count > (1U << 22)) {
        throw std::invalid_argument("Boozer result has too many modes");
    }
    const std::size_t points =
        checked_product(static_cast<std::size_t>(result.grid.ntheta),
                        static_cast<std::size_t>(result.grid.nzeta));
    const std::size_t real_count =
        checked_product(static_cast<std::size_t>(surfaces), points);
    const std::size_t spectral_count =
        checked_product(static_cast<std::size_t>(surfaces), mode_count);
    if (result.source_format_version < 1 || result.source_ns < 2 ||
        result.nfp < 1 || result.grid.first_surface < 1 ||
        result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        result.s.size() != static_cast<std::size_t>(surfaces) ||
        result.iota.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b2j00.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b.size() != real_count ||
        result.grid.sqrtg_b.size() != real_count ||
        result.spectrum.m.size() != mode_count ||
        result.spectrum.n.size() != mode_count ||
        result.spectrum.r.size() != spectral_count ||
        result.spectrum.z.size() != spectral_count ||
        result.spectrum.nu.size() != spectral_count) {
        throw std::invalid_argument("inconsistent Boozer result dimensions");
    }
    const auto finite_reals = [](std::span<const double> values) {
        return std::ranges::all_of(
            values, [](double value) { return std::isfinite(value); });
    };
    const auto finite_complex =
        [](std::span<const std::complex<double>> values) {
            return std::ranges::all_of(values, [](const auto value) {
                return std::isfinite(value.real()) &&
                       std::isfinite(value.imag());
            });
        };
    if (!finite_reals(result.s) || !finite_reals(result.iota) ||
        !finite_reals(result.grid.b) || !finite_reals(result.grid.sqrtg_b) ||
        !finite_reals(result.grid.b2j00) ||
        !finite_complex(result.spectrum.r) ||
        !finite_complex(result.spectrum.z) ||
        !finite_complex(result.spectrum.nu) ||
        !(result.resonance_tolerance > 0.0) ||
        !std::isfinite(result.resonance_tolerance)) {
        throw std::invalid_argument("Boozer result contains non-finite data");
    }
    std::size_t mode = 0;
    for (int n = -result.spectrum.nmax; n <= result.spectrum.nmax; ++n) {
        for (int m = -result.spectrum.mmax; m <= result.spectrum.mmax; ++m) {
            if (result.spectrum.m[mode] != m || result.spectrum.n[mode] != n) {
                throw std::invalid_argument(
                    "Boozer result has noncanonical mode ordering");
            }
            ++mode;
        }
    }
}

RealBoozerSpectrum make_real_spectrum(const BoozerResult& result) {
    validate_boozer_result(result);
    RealBoozerSpectrum output;
    output.mmax = result.spectrum.mmax;
    output.nmax = result.spectrum.nmax;
    const int real_modes = (output.mmax + 1) * (output.nmax + 1);
    const int signed_m_count = 2 * output.mmax + 1;
    const int signed_mode_count = signed_m_count * (2 * output.nmax + 1);
    const int surfaces = result.source_ns - result.grid.first_surface;
    output.m.reserve(static_cast<std::size_t>(real_modes));
    output.n.reserve(static_cast<std::size_t>(real_modes));
    for (int n = 0; n <= output.nmax; ++n) {
        for (int m = 0; m <= output.mmax; ++m) {
            output.m.push_back(m);
            output.n.push_back(n);
        }
    }
    const std::size_t coefficient_count = static_cast<std::size_t>(surfaces) *
                                          static_cast<std::size_t>(real_modes);
    output.rmncc.resize(coefficient_count);
    output.rmnss.resize(coefficient_count);
    output.zmnsc.resize(coefficient_count);
    output.zmncs.resize(coefficient_count);
    output.numnsc.resize(coefficient_count);
    output.numncs.resize(coefficient_count);

    const auto signed_mode = [&](int m, int n) {
        return (n + output.nmax) * signed_m_count + (m + output.mmax);
    };
    constexpr double FOUR_PI_SQUARED =
        4.0 * std::numbers::pi * std::numbers::pi;
    for (int surface = 0; surface < surfaces; ++surface) {
        for (int mode = 0; mode < real_modes; ++mode) {
            const int m = output.m[static_cast<std::size_t>(mode)];
            const int n = output.n[static_cast<std::size_t>(mode)];
            const std::size_t plus = static_cast<std::size_t>(
                surface * signed_mode_count + signed_mode(m, n));
            const std::size_t minus = static_cast<std::size_t>(
                surface * signed_mode_count + signed_mode(m, -n));
            const std::size_t target =
                static_cast<std::size_t>(surface * real_modes + mode);
            const double mode_weight = static_cast<double>(m == 0 ? 1 : 2) *
                                       static_cast<double>(n == 0 ? 1 : 2) /
                                       FOUR_PI_SQUARED;
            const auto r_plus = result.spectrum.r[plus];
            const auto r_minus = result.spectrum.r[minus];
            const auto z_plus = result.spectrum.z[plus];
            const auto z_minus = result.spectrum.z[minus];
            const auto nu_plus = result.spectrum.nu[plus];
            const auto nu_minus = result.spectrum.nu[minus];
            output.rmncc[target] =
                0.5 * mode_weight * (r_plus.real() + r_minus.real());
            output.rmnss[target] =
                0.5 * mode_weight * (r_minus.real() - r_plus.real());
            output.zmnsc[target] =
                -0.5 * mode_weight * (z_plus.imag() + z_minus.imag());
            output.zmncs[target] =
                0.5 * mode_weight * (z_minus.imag() - z_plus.imag());
            output.numnsc[target] =
                -0.5 * mode_weight * (nu_plus.imag() + nu_minus.imag());
            output.numncs[target] =
                0.5 * mode_weight * (nu_minus.imag() - nu_plus.imag());
        }
    }

    const auto restored = restore_complex_spectrum(output, result.source_ns,
                                                   result.grid.first_surface);
    const auto parity_matches = [](const auto& actual, const auto& expected) {
        for (std::size_t index = 0; index < actual.size(); ++index) {
            const double scale = std::max(1.0, std::abs(actual[index]));
            if (std::abs(actual[index] - expected[index]) > 1.0e-10 * scale) {
                return false;
            }
        }
        return true;
    };
    if (!parity_matches(result.spectrum.r, restored.r) ||
        !parity_matches(result.spectrum.z, restored.z) ||
        !parity_matches(result.spectrum.nu, restored.nu)) {
        throw std::invalid_argument(
            "Boozer spectra violate the real stellarator-symmetric parity "
            "contract");
    }
    return output;
}

MixedGridSpectrum restore_complex_spectrum(const RealBoozerSpectrum& spectrum,
                                           int source_ns,
                                           int first_surface) {
    const int surfaces = source_ns - first_surface;
    if (spectrum.mmax < 0 || spectrum.nmax < 0) {
        throw std::invalid_argument("invalid real Boozer mode bounds");
    }
    const std::size_t real_mode_count =
        checked_product(static_cast<std::size_t>(spectrum.mmax) + 1,
                        static_cast<std::size_t>(spectrum.nmax) + 1);
    const std::size_t signed_mode_count_size =
        checked_product(2 * static_cast<std::size_t>(spectrum.mmax) + 1,
                        2 * static_cast<std::size_t>(spectrum.nmax) + 1);
    if (real_mode_count > (1U << 20) || signed_mode_count_size > (1U << 22)) {
        throw std::invalid_argument("real Boozer spectrum has too many modes");
    }
    const int real_modes = static_cast<int>(real_mode_count);
    const int signed_m_count = 2 * spectrum.mmax + 1;
    const int signed_mode_count = static_cast<int>(signed_mode_count_size);
    const std::size_t expected = static_cast<std::size_t>(surfaces) *
                                 static_cast<std::size_t>(real_modes);
    if (surfaces < 1 ||
        spectrum.m.size() != static_cast<std::size_t>(real_modes) ||
        spectrum.n.size() != static_cast<std::size_t>(real_modes) ||
        spectrum.rmncc.size() != expected ||
        spectrum.rmnss.size() != expected ||
        spectrum.zmnsc.size() != expected ||
        spectrum.zmncs.size() != expected ||
        spectrum.numnsc.size() != expected ||
        spectrum.numncs.size() != expected) {
        throw std::invalid_argument("inconsistent real Boozer spectrum");
    }
    MixedGridSpectrum output;
    output.source_ns = source_ns;
    output.first_surface = first_surface;
    output.mmax = spectrum.mmax;
    output.nmax = spectrum.nmax;
    output.m.reserve(static_cast<std::size_t>(signed_mode_count));
    output.n.reserve(static_cast<std::size_t>(signed_mode_count));
    for (int n = -output.nmax; n <= output.nmax; ++n) {
        for (int m = -output.mmax; m <= output.mmax; ++m) {
            output.m.push_back(m);
            output.n.push_back(n);
        }
    }
    const std::size_t coefficient_count =
        static_cast<std::size_t>(surfaces) *
        static_cast<std::size_t>(signed_mode_count);
    output.r.assign(coefficient_count, {});
    output.z.assign(coefficient_count, {});
    output.nu.assign(coefficient_count, {});
    const auto signed_mode = [&](int m, int n) {
        return (n + output.nmax) * signed_m_count + (m + output.mmax);
    };
    constexpr double FOUR_PI_SQUARED =
        4.0 * std::numbers::pi * std::numbers::pi;
    for (int surface = 0; surface < surfaces; ++surface) {
        for (int mode = 0; mode < real_modes; ++mode) {
            const int m = spectrum.m[static_cast<std::size_t>(mode)];
            const int n = spectrum.n[static_cast<std::size_t>(mode)];
            const int expected_m = mode % (output.mmax + 1);
            const int expected_n = mode / (output.mmax + 1);
            if (m != expected_m || n != expected_n) {
                throw std::invalid_argument("noncanonical real Boozer modes");
            }
            const std::size_t source =
                static_cast<std::size_t>(surface * real_modes + mode);
            const double integral_scale =
                FOUR_PI_SQUARED / (static_cast<double>(m == 0 ? 1 : 2) *
                                   static_cast<double>(n == 0 ? 1 : 2));
            const std::complex<double> r_plus{
                integral_scale *
                    (spectrum.rmncc[source] - spectrum.rmnss[source]),
                0.0};
            const std::complex<double> r_minus{
                integral_scale *
                    (spectrum.rmncc[source] + spectrum.rmnss[source]),
                0.0};
            const std::complex<double> z_plus{
                0.0, -integral_scale *
                         (spectrum.zmnsc[source] + spectrum.zmncs[source])};
            const std::complex<double> z_minus{
                0.0, integral_scale *
                         (-spectrum.zmnsc[source] + spectrum.zmncs[source])};
            const std::complex<double> nu_plus{
                0.0, -integral_scale *
                         (spectrum.numnsc[source] + spectrum.numncs[source])};
            const std::complex<double> nu_minus{
                0.0, integral_scale *
                         (-spectrum.numnsc[source] + spectrum.numncs[source])};
            const auto target = [&](int signed_m, int signed_n) {
                return static_cast<std::size_t>(
                    surface * signed_mode_count +
                    signed_mode(signed_m, signed_n));
            };
            output.r[target(m, n)] = r_plus;
            output.r[target(m, -n)] = r_minus;
            output.r[target(-m, -n)] = std::conj(r_plus);
            output.r[target(-m, n)] = std::conj(r_minus);
            output.z[target(m, n)] = z_plus;
            output.z[target(m, -n)] = z_minus;
            output.z[target(-m, -n)] = std::conj(z_plus);
            output.z[target(-m, n)] = std::conj(z_minus);
            output.nu[target(m, n)] = nu_plus;
            output.nu[target(m, -n)] = nu_minus;
            output.nu[target(-m, -n)] = std::conj(nu_plus);
            output.nu[target(-m, n)] = std::conj(nu_minus);
        }
    }
    return output;
}

}  // namespace detail
}  // namespace magnetic_coordinate
