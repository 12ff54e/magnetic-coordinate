#ifndef MAGNETIC_COORDINATE_BOOZER_BINARY_HPP_
#define MAGNETIC_COORDINATE_BOOZER_BINARY_HPP_

#include "magnetic_coordinate/transform.hpp"

#include <filesystem>
#include <string>

namespace magnetic_coordinate {

struct BoozerFile {
    std::string source_path;
    BoozerResult result;
};

// Version 2 is little-endian and stores only real parity-family coefficients,
// matching cuMES's real spectral output style. Its coordinate metadata states
// that zeta is the unchanged source angle and that zeta_b = zeta + nu.
void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

BoozerFile read_boozer_binary(const std::filesystem::path& path);

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_BOOZER_BINARY_HPP_
