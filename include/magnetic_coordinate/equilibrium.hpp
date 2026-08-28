#ifndef MAGNETIC_COORDINATE_EQUILIBRIUM_HPP_
#define MAGNETIC_COORDINATE_EQUILIBRIUM_HPP_

#include <cstddef>
#include <span>
#include <vector>

namespace magnetic_coordinate {

// Non-owning view of the scientific fields in a cuMES equilibrium snapshot.
// Every array is [half_surface][zeta][theta], with theta contiguous.
struct NativeFieldView {
    int ns = 0;
    int ntheta = 0;
    int nzeta = 0;

    std::span<const double> sqrtg;
    std::span<const double> bsups;
    std::span<const double> bsupu;
    std::span<const double> bsupv;
    std::span<const double> bsubs;
    std::span<const double> bsubu;
    std::span<const double> bsubv;
};

// Fields and flux functions prepared on the ns-point integer radial grid.
// Real-space arrays retain the cuMES [surface][zeta][theta] layout.
struct IntegerGridEquilibrium {
    int ns = 0;
    int ntheta = 0;
    int nzeta = 0;

    std::vector<double> b2;
    std::vector<double> b;
    std::vector<double> sqrtg;

    // One value per integer radial surface.
    std::vector<double> toroidal_flux_derivative;
    std::vector<double> poloidal_flux_derivative;
    std::vector<double> iota;
};

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_EQUILIBRIUM_HPP_
