#ifndef MAGNETIC_COORDINATE_PEST_HPP_
#define MAGNETIC_COORDINATE_PEST_HPP_

#include "magnetic_coordinate/equilibrium.hpp"
#include "magnetic_coordinate/spectral.hpp"

#include <vector>

namespace magnetic_coordinate {

// Regular theta_p/original-zeta grid for non-axis surfaces. Surface index zero
// in these arrays corresponds to source integer surface first_surface.
struct PestGrid {
    int source_ns = 0;
    int first_surface = 1;
    int ntheta = 0;
    int nzeta = 0;

    std::vector<double> source_theta;
    std::vector<double> b2;
    std::vector<double> sqrtg;
};

PestGrid remap_to_pest(const NativeAngularGeometry& geometry,
                       const IntegerGridEquilibrium& integer_fields);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_PEST_HPP_
