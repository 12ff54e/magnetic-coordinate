#ifndef MAGNETIC_COORDINATE_BOOZER_FFT_HPP_
#define MAGNETIC_COORDINATE_BOOZER_FFT_HPP_

#include "magnetic_coordinate/pest.hpp"

#include <span>
#include <vector>

namespace magnetic_coordinate {

struct BoozerShift {
    int source_ns = 0;
    int first_surface = 1;
    int ntheta = 0;
    int nzeta = 0;

    // nu(theta_p,zeta), surface-major with theta contiguous.
    std::vector<double> nu;

    // Continuous Fourier-integral zero mode of B^2*sqrt(g_p):
    // b2j00 = 4*pi^2*angular_mean(B^2*sqrt(g_p)).
    std::vector<double> b2j00;
};

// Batched GPU transform and spectral magnetic differential-equation solve.
// `iota` is on all source integer surfaces; surface zero is skipped with the
// PEST axis. A nonzero resonant numerator is reported as a domain error.
BoozerShift solve_boozer_shift_gpu(const PestGrid& pest,
                                   std::span<const double> iota,
                                   int nfp,
                                   double resonance_tolerance = 1.0e-12);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_BOOZER_FFT_HPP_
