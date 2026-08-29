#include "magnetic_coordinate/boozer_binary.hpp"

#include "boozer_output_internal.hpp"

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace magnetic_coordinate {
namespace {

constexpr std::array<char, 8> MAGIC{'M', 'C', 'B', 'O', 'O', 'Z', '0', '2'};
constexpr std::int32_t VERSION = 2;

std::size_t checked_product(std::size_t left, std::size_t right) {
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left) {
        throw std::runtime_error("Boozer binary dimension overflow");
    }
    return left * right;
}

class Writer {
   public:
    explicit Writer(const std::filesystem::path& path)
        : output_(path, std::ios::binary | std::ios::trunc) {
        if (!output_) throw std::runtime_error("cannot open Boozer output");
    }
    void bytes(const void* data, std::size_t count) {
        output_.write(static_cast<const char*>(data),
                      static_cast<std::streamsize>(count));
        if (!output_) throw std::runtime_error("could not write Boozer output");
    }
    template <typename UInt>
    void uint(UInt value) {
        static_assert(std::is_unsigned_v<UInt>);
        std::array<std::uint8_t, sizeof(UInt)> data{};
        for (std::size_t i = 0; i < data.size(); ++i) {
            data[i] = static_cast<std::uint8_t>(value >> (8U * i));
        }
        bytes(data.data(), data.size());
    }
    void i32(std::int32_t value) { uint(std::bit_cast<std::uint32_t>(value)); }
    void f64(double value) { uint(std::bit_cast<std::uint64_t>(value)); }
    void string(std::string_view value) {
        if (value.size() > static_cast<std::size_t>(
                               std::numeric_limits<std::int32_t>::max())) {
            throw std::runtime_error("Boozer output string is too large");
        }
        i32(static_cast<std::int32_t>(value.size()));
        bytes(value.data(), value.size());
    }
    void doubles(std::span<const double> values) {
        for (double value : values) f64(value);
    }
    void integers(std::span<const int> values) {
        for (int value : values) i32(value);
    }

   private:
    std::ofstream output_;
};

class Reader {
   public:
    explicit Reader(const std::filesystem::path& path)
        : input_(path, std::ios::binary) {
        if (!input_) throw std::runtime_error("cannot open Boozer input");
    }
    void bytes(void* data, std::size_t count) {
        input_.read(static_cast<char*>(data),
                    static_cast<std::streamsize>(count));
        if (!input_) throw std::runtime_error("truncated Boozer input");
    }
    template <typename UInt>
    UInt uint() {
        static_assert(std::is_unsigned_v<UInt>);
        std::array<std::uint8_t, sizeof(UInt)> data{};
        bytes(data.data(), data.size());
        UInt value = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            value |= static_cast<UInt>(data[i]) << (8U * i);
        }
        return value;
    }
    std::int32_t i32() {
        return std::bit_cast<std::int32_t>(uint<std::uint32_t>());
    }
    double f64() { return std::bit_cast<double>(uint<std::uint64_t>()); }
    std::string string() {
        const auto count = i32();
        if (count < 0 || count > (1 << 24)) {
            throw std::runtime_error("invalid Boozer string length");
        }
        std::string value(static_cast<std::size_t>(count), '\0');
        bytes(value.data(), value.size());
        return value;
    }
    std::vector<double> doubles(std::size_t count) {
        std::vector<double> values(count);
        for (double& value : values) value = f64();
        return values;
    }
    std::vector<int> integers(std::size_t count) {
        std::vector<int> values(count);
        for (int& value : values) value = i32();
        return values;
    }

   private:
    std::ifstream input_;
};

}  // namespace

void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    const detail::RealBoozerSpectrum spectrum =
        detail::make_real_spectrum(result);
    Writer writer(path);
    writer.bytes(MAGIC.data(), MAGIC.size());
    writer.i32(VERSION);
    writer.string(detail::COORDINATE_CONVENTION);
    writer.string(detail::FOURIER_CONVENTION);
    writer.string(source_path.string());
    writer.i32(result.source_format_version);
    writer.i32(result.source_ns);
    writer.i32(result.source_ntheta);
    writer.i32(result.source_nzeta);
    writer.i32(result.source_mpol);
    writer.i32(result.source_ntor);
    writer.i32(result.nfp);
    writer.i32(result.grid.first_surface);
    writer.i32(result.grid.ntheta);
    writer.i32(result.grid.nzeta);
    writer.i32(result.spectrum.mmax);
    writer.i32(result.spectrum.nmax);
    writer.i32(static_cast<std::int32_t>(result.radial_order));
    writer.f64(result.resonance_tolerance);
    writer.doubles(result.s);
    writer.doubles(result.iota);
    writer.integers(spectrum.m);
    writer.integers(spectrum.n);
    writer.doubles(spectrum.rmncc);
    writer.doubles(spectrum.rmnss);
    writer.doubles(spectrum.zmnsc);
    writer.doubles(spectrum.zmncs);
    writer.doubles(spectrum.numnsc);
    writer.doubles(spectrum.numncs);
    writer.doubles(result.grid.b);
    writer.doubles(result.grid.sqrtg_b);
    writer.doubles(result.grid.b2j00);
}

BoozerFile read_boozer_binary(const std::filesystem::path& path) {
    Reader reader(path);
    std::array<char, 8> magic{};
    reader.bytes(magic.data(), magic.size());
    if (magic != MAGIC || reader.i32() != VERSION) {
        throw std::runtime_error("unsupported Boozer binary format");
    }
    if (reader.string() != detail::COORDINATE_CONVENTION ||
        reader.string() != detail::FOURIER_CONVENTION) {
        throw std::runtime_error("unsupported Boozer coordinate convention");
    }
    BoozerFile file;
    file.source_path = reader.string();
    auto& result = file.result;
    result.source_format_version = reader.i32();
    result.source_ns = reader.i32();
    result.source_ntheta = reader.i32();
    result.source_nzeta = reader.i32();
    result.source_mpol = reader.i32();
    result.source_ntor = reader.i32();
    result.nfp = reader.i32();
    result.grid.first_surface = reader.i32();
    result.grid.ntheta = reader.i32();
    result.grid.nzeta = reader.i32();
    result.spectrum.mmax = reader.i32();
    result.spectrum.nmax = reader.i32();
    result.radial_order = static_cast<RadialInterpolationOrder>(reader.i32());
    result.resonance_tolerance = reader.f64();
    result.grid.source_ns = result.source_ns;
    const int surfaces = result.source_ns - result.grid.first_surface;
    if (result.spectrum.mmax < 0 || result.spectrum.nmax < 0 ||
        result.spectrum.mmax > (1 << 15) || result.spectrum.nmax > (1 << 15)) {
        throw std::runtime_error("invalid Boozer binary mode bounds");
    }
    const std::size_t mode_count =
        checked_product(static_cast<std::size_t>(result.spectrum.mmax) + 1,
                        static_cast<std::size_t>(result.spectrum.nmax) + 1);
    if (surfaces < 1 || result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        mode_count < 1 || surfaces > (1 << 20) || mode_count > (1 << 20)) {
        throw std::runtime_error("invalid Boozer binary dimensions");
    }
    const std::size_t real_count = checked_product(
        static_cast<std::size_t>(surfaces),
        checked_product(static_cast<std::size_t>(result.grid.ntheta),
                        static_cast<std::size_t>(result.grid.nzeta)));
    const std::size_t coefficient_count =
        checked_product(static_cast<std::size_t>(surfaces), mode_count);
    result.s = reader.doubles(static_cast<std::size_t>(surfaces));
    result.iota = reader.doubles(static_cast<std::size_t>(surfaces));
    detail::RealBoozerSpectrum spectrum;
    spectrum.mmax = result.spectrum.mmax;
    spectrum.nmax = result.spectrum.nmax;
    spectrum.m = reader.integers(mode_count);
    spectrum.n = reader.integers(mode_count);
    spectrum.rmncc = reader.doubles(coefficient_count);
    spectrum.rmnss = reader.doubles(coefficient_count);
    spectrum.zmnsc = reader.doubles(coefficient_count);
    spectrum.zmncs = reader.doubles(coefficient_count);
    spectrum.numnsc = reader.doubles(coefficient_count);
    spectrum.numncs = reader.doubles(coefficient_count);
    result.spectrum = detail::restore_complex_spectrum(
        spectrum, result.source_ns, result.grid.first_surface);
    result.grid.b = reader.doubles(real_count);
    result.grid.sqrtg_b = reader.doubles(real_count);
    result.grid.b2j00 = reader.doubles(static_cast<std::size_t>(surfaces));
    detail::validate_boozer_result(result);
    return file;
}

}  // namespace magnetic_coordinate
