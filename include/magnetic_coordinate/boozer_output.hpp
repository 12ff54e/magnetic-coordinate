#ifndef MAGNETIC_COORDINATE_BOOZER_OUTPUT_HPP_
#define MAGNETIC_COORDINATE_BOOZER_OUTPUT_HPP_

#include "magnetic_coordinate/transform.hpp"

#include <cstdint>
#include <filesystem>
#include <string_view>

namespace magnetic_coordinate {

enum class BoozerOutputFormat : std::uint8_t {
    BINARY = 0,
    NETCDF = 1,
    HDF5 = 2,
};

struct BoozerOutputSpec {
    BoozerOutputFormat format = BoozerOutputFormat::BINARY;
    std::filesystem::path path;
};

// Resolve .bin, .nc, .h5, or .hdf5 case-insensitively. Unknown and missing
// suffixes are rejected instead of receiving misleading binary contents.
BoozerOutputSpec resolve_boozer_output_spec(const std::filesystem::path& path);

bool boozer_output_format_available(BoozerOutputFormat format);
std::string_view boozer_output_suffix(BoozerOutputFormat format);

// Write the selected v3 container. The result remains format-neutral and can
// instead be retained in memory by a future solver-library caller.
void write_boozer_output(const BoozerOutputSpec& spec,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path);

inline void write_boozer_output(const std::filesystem::path& path,
                                const BoozerResult& result,
                                const std::filesystem::path& source_path) {
    write_boozer_output(resolve_boozer_output_spec(path), result, source_path);
}

}  // namespace magnetic_coordinate

#endif  // MAGNETIC_COORDINATE_BOOZER_OUTPUT_HPP_
