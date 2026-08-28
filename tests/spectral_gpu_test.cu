#include "magnetic_coordinate/spectral.hpp"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

void test_gpu_matches_cpu() {
    magnetic_coordinate::CumesEquilibrium equilibrium;
    equilibrium.ns = 5;
    equilibrium.mpol = 4;
    equilibrium.ntor = 2;
    equilibrium.nfp = 5;
    equilibrium.mnmax = equilibrium.mpol * (equilibrium.ntor + 1);
    equilibrium.ntheta = 12;
    equilibrium.nzeta = 8;
    equilibrium.phiedge = -1.0;
    equilibrium.aphi = {1.0};
    const std::size_t family_size =
        static_cast<std::size_t>(equilibrium.ns * equilibrium.mnmax);
    for (std::size_t family = 0;
         family < magnetic_coordinate::CumesEquilibrium::SPECTRAL_FAMILY_COUNT;
         ++family) {
        equilibrium.families[family].resize(family_size);
        for (std::size_t coefficient = 0; coefficient < family_size;
             ++coefficient) {
            equilibrium.families[family][coefficient] =
                0.01 * static_cast<double>(family + 1) *
                std::sin(0.3 * static_cast<double>(coefficient + 1));
        }
    }
    const auto normalization =
        magnetic_coordinate::reconstruct_flux_normalization(equilibrium);
    const auto cpu = magnetic_coordinate::synthesize_native_geometry(
        equilibrium, normalization);
    const auto gpu = magnetic_coordinate::synthesize_native_geometry_gpu(
        equilibrium, normalization);
    const auto compare = [](const std::vector<double>& left,
                            const std::vector<double>& right) {
        if (left.size() != right.size()) {
            throw std::runtime_error("GPU synthesis extent mismatch");
        }
        double largest = 0.0;
        for (std::size_t index = 0; index < left.size(); ++index) {
            largest = std::max(largest, std::abs(left[index] - right[index]));
        }
        if (largest > 2.0e-14) {
            throw std::runtime_error("GPU synthesis differs from CPU oracle");
        }
    };
    compare(cpu.r, gpu.r);
    compare(cpu.z, gpu.z);
    compare(cpu.lambda, gpu.lambda);
    compare(cpu.lambda_theta, gpu.lambda_theta);
    compare(cpu.lambda_zeta, gpu.lambda_zeta);
}

}  // namespace

int main() {
    int devices = 0;
    const cudaError_t status = cudaGetDeviceCount(&devices);
    if (status == cudaErrorNoDevice || status == cudaErrorInsufficientDriver ||
        devices == 0) {
        std::cout << "spectral_gpu_test: SKIP (no CUDA device)\n";
        return 0;
    }
    if (status != cudaSuccess) {
        std::cerr << "spectral_gpu_test: " << cudaGetErrorString(status) << '\n';
        return 1;
    }
    try {
        test_gpu_matches_cpu();
    } catch (const std::exception& error) {
        std::cerr << "spectral_gpu_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "spectral_gpu_test: PASS\n";
    return 0;
}
