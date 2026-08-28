#ifndef MAGNETIC_COORDINATE_MIXED_SPECTRUM_HPP_
#define MAGNETIC_COORDINATE_MIXED_SPECTRUM_HPP_

#include "magnetic_coordinate/mixed_grid.hpp"

#include <complex>
#include <vector>

namespace magnetic_coordinate {

struct MixedGridSpectrum {
    int source_ns = 0;
    int first_surface = 1;
    int mmax = 0;
    int nmax = 0;

    // Mode ordering is n-major, then m from -mmax through mmax.
    std::vector<int> m;
    std::vector<int> n;
    std::vector<std::complex<double>> r;
    std::vector<std::complex<double>> z;
    std::vector<std::complex<double>> nu;
};

// Batched GPU forward analysis. Coefficients are continuous angular-integral
// coefficients: raw unnormalized FFT values times 4*pi^2/(ntheta*nzeta).
MixedGridSpectrum analyze_mixed_grid_gpu(const BoozerMixedGrid& grid,
                                         int mmax,
                                         int nmax);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_MIXED_SPECTRUM_HPP_
