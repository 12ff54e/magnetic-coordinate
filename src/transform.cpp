#include "magnetic_coordinate/transform.hpp"

#include "magnetic_coordinate/angular_resample.hpp"
#include "magnetic_coordinate/boozer_fft.hpp"
#include "magnetic_coordinate/pest.hpp"
#include "magnetic_coordinate/spectral.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace magnetic_coordinate {

BoozerResult transform_to_boozer(const CumesEquilibrium& equilibrium,
                                 const TransformSettings& settings) {
    if (settings.output_ntheta < 0 || settings.output_nzeta < 0 ||
        settings.mmax < -1 ||
        settings.nmax < -1) {
        throw std::invalid_argument("invalid Boozer transform settings");
    }
    const int output_ntheta = settings.output_ntheta == 0
                                  ? equilibrium.ntheta
                                  : settings.output_ntheta;
    const int output_nzeta = settings.output_nzeta == 0
                                 ? equilibrium.nzeta
                                 : settings.output_nzeta;
    const int available_m = (output_ntheta - 1) / 2;
    const int available_n = (output_nzeta - 1) / 2;
    const int mmax = settings.mmax < 0
                         ? std::min(equilibrium.mpol - 1, available_m)
                         : settings.mmax;
    const int nmax = settings.nmax < 0
                         ? std::min(equilibrium.ntor, available_n)
                         : settings.nmax;

    const auto source_integer_fields = prepare_integer_grid(
        equilibrium.native_field_view(), settings.radial_order);
    const auto normalization = reconstruct_flux_normalization(equilibrium);
    const auto source_geometry =
        synthesize_native_geometry_gpu(equilibrium, normalization);
    const auto angular = resample_toroidal_grid(
        source_geometry, source_integer_fields, output_nzeta);
    const auto pest = remap_to_pest(angular.geometry, angular.fields);
    const auto shift = solve_boozer_shift_gpu(
        pest, angular.fields.iota, equilibrium.nfp,
        settings.resonance_tolerance);
    auto grid = remap_to_boozer_mixed_grid(
        angular.geometry, pest, shift, angular.fields.iota, output_ntheta);
    auto spectrum = analyze_mixed_grid_gpu(grid, mmax, nmax);

    BoozerResult result;
    result.source_format_version = equilibrium.format_version;
    result.source_ns = equilibrium.ns;
    result.source_ntheta = equilibrium.ntheta;
    result.source_nzeta = equilibrium.nzeta;
    result.source_mpol = equilibrium.mpol;
    result.source_ntor = equilibrium.ntor;
    result.nfp = equilibrium.nfp;
    result.radial_order = settings.radial_order;
    result.resonance_tolerance = settings.resonance_tolerance;
    const int surfaces = equilibrium.ns - grid.first_surface;
    result.s.resize(static_cast<std::size_t>(surfaces));
    result.iota.resize(static_cast<std::size_t>(surfaces));
    for (int surface = 0; surface < surfaces; ++surface) {
        const int source_surface = surface + grid.first_surface;
        result.s[static_cast<std::size_t>(surface)] =
            static_cast<double>(source_surface) /
            static_cast<double>(equilibrium.ns - 1);
        result.iota[static_cast<std::size_t>(surface)] =
            angular.fields.iota[static_cast<std::size_t>(source_surface)];
    }
    result.grid = std::move(grid);
    result.spectrum = std::move(spectrum);
    return result;
}

BoozerResult transform_cumes_file(const std::filesystem::path& input,
                                  const TransformSettings& settings) {
    return transform_to_boozer(read_cumes_binary(input), settings);
}

}  // namespace magnetic_coordinate
