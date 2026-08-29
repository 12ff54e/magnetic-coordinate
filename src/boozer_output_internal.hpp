#ifndef MAGNETIC_COORDINATE_SRC_BOOZER_OUTPUT_INTERNAL_HPP_
#define MAGNETIC_COORDINATE_SRC_BOOZER_OUTPUT_INTERNAL_HPP_

#include "magnetic_coordinate/transform.hpp"

#include <filesystem>
#include <string_view>
#include <vector>

namespace magnetic_coordinate::detail {

inline constexpr std::string_view BOOZER_SCHEMA =
    "magnetic-coordinate-boozer-v2";
inline constexpr std::string_view COORDINATE_CONVENTION =
    "mixed-grid-v1: theta_b uniform; zeta is the unchanged source toroidal "
    "angle; zeta_b=zeta+nu";
inline constexpr std::string_view FOURIER_CONVENTION =
    "real-parity-v2: f=sum[cc*cos(m*theta_b)*cos(n*zeta) + "
    "ss*sin(m*theta_b)*sin(n*zeta)] for even fields and "
    "f=sum[sc*sin(m*theta_b)*cos(n*zeta) + "
    "cs*cos(m*theta_b)*sin(n*zeta)] for odd fields; m,n are nonnegative "
    "and n is a field-period mode";

struct RealBoozerSpectrum {
    int mmax = 0;
    int nmax = 0;
    std::vector<int> m;
    std::vector<int> n;
    std::vector<double> rmncc;
    std::vector<double> rmnss;
    std::vector<double> zmnsc;
    std::vector<double> zmncs;
    std::vector<double> numnsc;
    std::vector<double> numncs;
};

void validate_boozer_result(const BoozerResult& result);
RealBoozerSpectrum make_real_spectrum(const BoozerResult& result);
MixedGridSpectrum restore_complex_spectrum(const RealBoozerSpectrum& spectrum,
                                           int source_ns,
                                           int first_surface);

#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
void write_boozer_netcdf(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);
#endif

#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
void write_boozer_hdf5(const std::filesystem::path& path,
                       const BoozerResult& result,
                       const std::filesystem::path& source_path);
#endif

}  // namespace magnetic_coordinate::detail

#endif  // MAGNETIC_COORDINATE_SRC_BOOZER_OUTPUT_INTERNAL_HPP_
