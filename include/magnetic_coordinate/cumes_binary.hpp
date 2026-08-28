#ifndef MAGNETIC_COORDINATE_CUMES_BINARY_HPP_
#define MAGNETIC_COORDINATE_CUMES_BINARY_HPP_

#include "magnetic_coordinate/equilibrium.hpp"

#include <array>
#include <cstdint>
#include <filesystem>
#include <vector>

namespace magnetic_coordinate {

struct CumesEquilibrium {
    static constexpr std::size_t SPECTRAL_FAMILY_COUNT = 6;
    static constexpr std::size_t HALF_FIELD_COUNT = 7;
    static constexpr std::size_t FULL_FIELD_COUNT = 6;

    enum SpectralFamily : std::uint8_t {
        RMNCC = 0,
        ZMNSC = 1,
        LMNSC = 2,
        RMNSS = 3,
        ZMNCS = 4,
        LMNCS = 5,
    };

    enum HalfField : std::uint8_t {
        SQRTG = 0,
        BSUPS = 1,
        BSUPU = 2,
        BSUPV = 3,
        BSUBS = 4,
        BSUBU = 5,
        BSUBV = 6,
    };

    int format_version = 0;
    int ns = 0;
    int mnmax = 0;
    int mpol = 0;
    int ntor = 0;
    int nfp = 0;
    int ntheta = 0;
    int nzeta = 0;
    int ncurr = 0;
    double phiedge = 0.0;
    std::vector<double> aphi;

    std::array<std::vector<double>, SPECTRAL_FAMILY_COUNT> families;
    std::array<std::vector<double>, HALF_FIELD_COUNT> half_fields;
    std::array<std::vector<double>, FULL_FIELD_COUNT> full_fields;

    NativeFieldView native_field_view() const;
};

// Read the current cuMES schema-v1 binary result. The magnetic-coordinate
// transform requires the scientific fields introduced by on-disk version 8,
// so older state-only versions are rejected deliberately.
CumesEquilibrium read_cumes_binary(const std::filesystem::path& path);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_CUMES_BINARY_HPP_
