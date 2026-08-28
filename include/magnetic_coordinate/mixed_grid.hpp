#ifndef MAGNETIC_COORDINATE_MIXED_GRID_HPP_
#define MAGNETIC_COORDINATE_MIXED_GRID_HPP_

#include "magnetic_coordinate/boozer_fft.hpp"
#include "magnetic_coordinate/spectral.hpp"

#include <span>
#include <vector>

namespace magnetic_coordinate {

// Fields on a grid uniform in Boozer poloidal angle while retaining the
// source toroidal angle. Surface zero corresponds to first_surface in the
// source equilibrium; theta is contiguous in every array.
struct BoozerMixedGrid {
    int source_ns = 0;
    int first_surface = 1;
    int ntheta = 0;
    int nzeta = 0;

    std::vector<double> source_theta_p;
    std::vector<double> source_theta;
    std::vector<double> r;
    std::vector<double> z;
    std::vector<double> nu;
    std::vector<double> b;
    std::vector<double> sqrtg_b;
    std::vector<double> b2j00;
};

// The poloidal resolution is independent of the source grid. The toroidal
// coordinate deliberately remains the unchanged source zeta grid.
BoozerMixedGrid remap_to_boozer_mixed_grid(
    const NativeAngularGeometry& geometry,
    const PestGrid& pest,
    const BoozerShift& shift,
    std::span<const double> iota,
    int output_ntheta);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_MIXED_GRID_HPP_
