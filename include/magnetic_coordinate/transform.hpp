#ifndef MAGNETIC_COORDINATE_TRANSFORM_HPP_
#define MAGNETIC_COORDINATE_TRANSFORM_HPP_

#include "magnetic_coordinate/cumes_binary.hpp"
#include "magnetic_coordinate/mixed_spectrum.hpp"
#include "magnetic_coordinate/radial.hpp"

#include <filesystem>
#include <vector>

namespace magnetic_coordinate {

struct TransformSettings {
    // Zero selects the source angular grid size. Negative truncations select
    // the source truncation, clipped below the output Nyquist modes.
    int output_ntheta = 0;
    int mmax = -1;
    int nmax = -1;
    RadialInterpolationOrder radial_order =
        RadialInterpolationOrder::FOUR_POINT;
    double resonance_tolerance = 1.0e-12;
};

struct BoozerResult {
    int source_format_version = 0;
    int source_ns = 0;
    int source_ntheta = 0;
    int source_nzeta = 0;
    int source_mpol = 0;
    int source_ntor = 0;
    int nfp = 0;
    RadialInterpolationOrder radial_order =
        RadialInterpolationOrder::FOUR_POINT;
    double resonance_tolerance = 1.0e-12;

    // Normalized toroidal-flux coordinate and iota on exported non-axis
    // surfaces.
    std::vector<double> s;
    std::vector<double> iota;
    BoozerMixedGrid grid;
    MixedGridSpectrum spectrum;
};

BoozerResult transform_to_boozer(const CumesEquilibrium& equilibrium,
                                 const TransformSettings& settings = {});

BoozerResult transform_cumes_file(
    const std::filesystem::path& input,
    const TransformSettings& settings = {});

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_TRANSFORM_HPP_
