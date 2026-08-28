#ifndef MAGNETIC_COORDINATE_SPECTRAL_HPP_
#define MAGNETIC_COORDINATE_SPECTRAL_HPP_

#include "magnetic_coordinate/cumes_binary.hpp"

#include <vector>

namespace magnetic_coordinate {

struct FluxNormalization {
    double maximum_toroidal_flux = 0.0;
    double lambda_scale = 0.0;
    std::vector<double> phip_full;
    std::vector<double> phip_half;
    std::vector<double> lambda_multiplier;
};

// Reproduce the cuMES toroidal-flux profile and lambda normalization. The
// physical PEST displacement is lambda_multiplier[s] times the stored lambda
// spectral state on each integer surface.
FluxNormalization reconstruct_flux_normalization(
    const CumesEquilibriumView& equilibrium);

inline FluxNormalization reconstruct_flux_normalization(
    const CumesEquilibrium& equilibrium) {
    return reconstruct_flux_normalization(equilibrium.view());
}

struct NativeAngularGeometry {
    int ns = 0;
    int ntheta = 0;
    int nzeta = 0;

    std::vector<double> r;
    std::vector<double> z;
    std::vector<double> lambda;
    std::vector<double> lambda_theta;
    std::vector<double> lambda_zeta;
};

// CPU reference synthesis of the folded cuMES spectral families. The zeta
// derivative is with respect to the physical toroidal angle and therefore
// includes nfp. GPU synthesis must reproduce this reference.
NativeAngularGeometry synthesize_native_geometry(
    const CumesEquilibriumView& equilibrium,
    const FluxNormalization& normalization);

inline NativeAngularGeometry synthesize_native_geometry(
    const CumesEquilibrium& equilibrium,
    const FluxNormalization& normalization) {
    return synthesize_native_geometry(equilibrium.view(), normalization);
}

// Production synthesis across every surface and angular point on the GPU.
// The CPU implementation above remains the independent reference oracle.
NativeAngularGeometry synthesize_native_geometry_gpu(
    const CumesEquilibriumView& equilibrium,
    const FluxNormalization& normalization);

inline NativeAngularGeometry synthesize_native_geometry_gpu(
    const CumesEquilibrium& equilibrium,
    const FluxNormalization& normalization) {
    return synthesize_native_geometry_gpu(equilibrium.view(), normalization);
}

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_SPECTRAL_HPP_
