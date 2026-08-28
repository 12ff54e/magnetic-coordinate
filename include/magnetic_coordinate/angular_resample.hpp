#ifndef MAGNETIC_COORDINATE_ANGULAR_RESAMPLE_HPP_
#define MAGNETIC_COORDINATE_ANGULAR_RESAMPLE_HPP_

#include "magnetic_coordinate/equilibrium.hpp"
#include "magnetic_coordinate/spectral.hpp"

namespace magnetic_coordinate {

struct ResampledAngularData {
    NativeAngularGeometry geometry;
    IntegerGridEquilibrium fields;
};

// Periodically resample the unchanged field-period zeta coordinate. The
// poloidal grid and radial surfaces are not modified. Axisymmetric data are
// replicated exactly.
ResampledAngularData resample_toroidal_grid(
    const NativeAngularGeometry& geometry,
    const IntegerGridEquilibrium& fields,
    int output_nzeta);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_ANGULAR_RESAMPLE_HPP_
