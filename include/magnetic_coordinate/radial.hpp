#ifndef MAGNETIC_COORDINATE_RADIAL_HPP_
#define MAGNETIC_COORDINATE_RADIAL_HPP_

#include "magnetic_coordinate/equilibrium.hpp"

#include <cstddef>
#include <span>
#include <vector>

namespace magnetic_coordinate {

enum class RadialInterpolationOrder {
    TWO_POINT = 2,
    FOUR_POINT = 4,
};

// Interpolate values from s=(j+1/2)/(ns-1), j=0..ns-2, to
// s=j/(ns-1), j=0..ns-1. Endpoint values use the same Lagrange polynomial as
// the interior and are therefore explicit extrapolations.
std::vector<double> interpolate_half_to_full(
    std::span<const double> half_values,
    int ns,
    std::size_t values_per_surface,
    RadialInterpolationOrder order);

// Validate native cuMES fields, compute B^2=B^i B_i on the half grid,
// unstagger B^2 and sqrt(g), and recover the converged flux functions from
// angular averages of sqrt(g) B^theta and sqrt(g) B^zeta.
IntegerGridEquilibrium prepare_integer_grid(
    const NativeFieldView& native,
    RadialInterpolationOrder order = RadialInterpolationOrder::FOUR_POINT);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_RADIAL_HPP_
