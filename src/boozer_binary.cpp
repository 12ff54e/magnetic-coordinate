#include "magnetic_coordinate/boozer_binary.hpp"

#include <array>
#include <bit>
#include <complex>
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

constexpr std::array<char, 8> MAGIC{'M', 'C', 'B', 'O', 'O', 'Z', '0', '1'};
constexpr std::int32_t VERSION = 1;
constexpr std::string_view COORDINATE_CONVENTION =
    "mixed-grid-v1: theta_b uniform; zeta is the unchanged source toroidal "
    "angle; zeta_b=zeta+nu";
constexpr std::string_view FOURIER_CONVENTION =
    "integral-v1: integral f*exp(-i(m*theta_b+n*zeta)) dtheta_b dzeta; "
    "n is a field-period mode and physical derivatives multiply n by nfp";

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
        if (value.size() >
            static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
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
    void complex(std::span<const std::complex<double>> values) {
        for (const auto value : values) {
            f64(value.real());
            f64(value.imag());
        }
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
        input_.read(static_cast<char*>(data), static_cast<std::streamsize>(count));
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
    std::int32_t i32() { return std::bit_cast<std::int32_t>(uint<std::uint32_t>()); }
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
    std::vector<std::complex<double>> complex(std::size_t count) {
        std::vector<std::complex<double>> values(count);
        for (auto& value : values) value = {f64(), f64()};
        return values;
    }

   private:
    std::ifstream input_;
};

void validate(const BoozerResult& result) {
    const int surfaces = result.source_ns - result.grid.first_surface;
    const int mode_count =
        (2 * result.spectrum.mmax + 1) * (2 * result.spectrum.nmax + 1);
    const std::size_t real_count = checked_product(
        static_cast<std::size_t>(surfaces),
        static_cast<std::size_t>(result.grid.ntheta * result.grid.nzeta));
    const std::size_t spectral_count = checked_product(
        static_cast<std::size_t>(surfaces),
        static_cast<std::size_t>(mode_count));
    if (result.source_format_version < 1 || result.source_ns < 2 ||
        result.nfp < 1 || result.grid.first_surface < 1 ||
        result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        result.s.size() != static_cast<std::size_t>(surfaces) ||
        result.iota.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b2j00.size() != static_cast<std::size_t>(surfaces) ||
        result.grid.b.size() != real_count ||
        result.grid.sqrtg_b.size() != real_count ||
        result.spectrum.m.size() != static_cast<std::size_t>(mode_count) ||
        result.spectrum.n.size() != static_cast<std::size_t>(mode_count) ||
        result.spectrum.r.size() != spectral_count ||
        result.spectrum.z.size() != spectral_count ||
        result.spectrum.nu.size() != spectral_count) {
        throw std::invalid_argument("inconsistent Boozer result dimensions");
    }
}

}  // namespace

void write_boozer_binary(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    validate(result);
    Writer writer(path);
    writer.bytes(MAGIC.data(), MAGIC.size());
    writer.i32(VERSION);
    writer.string(COORDINATE_CONVENTION);
    writer.string(FOURIER_CONVENTION);
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
    writer.integers(result.spectrum.m);
    writer.integers(result.spectrum.n);
    writer.complex(result.spectrum.r);
    writer.complex(result.spectrum.z);
    writer.complex(result.spectrum.nu);
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
    if (reader.string() != COORDINATE_CONVENTION ||
        reader.string() != FOURIER_CONVENTION) {
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
    result.radial_order =
        static_cast<RadialInterpolationOrder>(reader.i32());
    result.resonance_tolerance = reader.f64();
    result.grid.source_ns = result.source_ns;
    result.spectrum.source_ns = result.source_ns;
    result.spectrum.first_surface = result.grid.first_surface;
    const int surfaces = result.source_ns - result.grid.first_surface;
    const int mode_count =
        (2 * result.spectrum.mmax + 1) * (2 * result.spectrum.nmax + 1);
    if (surfaces < 1 || result.grid.ntheta < 1 || result.grid.nzeta < 1 ||
        mode_count < 1 || surfaces > (1 << 20) || mode_count > (1 << 20)) {
        throw std::runtime_error("invalid Boozer binary dimensions");
    }
    const std::size_t real_count = checked_product(
        static_cast<std::size_t>(surfaces),
        checked_product(static_cast<std::size_t>(result.grid.ntheta),
                        static_cast<std::size_t>(result.grid.nzeta)));
    const std::size_t spectral_count = checked_product(
        static_cast<std::size_t>(surfaces),
        static_cast<std::size_t>(mode_count));
    result.s = reader.doubles(static_cast<std::size_t>(surfaces));
    result.iota = reader.doubles(static_cast<std::size_t>(surfaces));
    result.spectrum.m = reader.integers(static_cast<std::size_t>(mode_count));
    result.spectrum.n = reader.integers(static_cast<std::size_t>(mode_count));
    result.spectrum.r = reader.complex(spectral_count);
    result.spectrum.z = reader.complex(spectral_count);
    result.spectrum.nu = reader.complex(spectral_count);
    result.grid.b = reader.doubles(real_count);
    result.grid.sqrtg_b = reader.doubles(real_count);
    result.grid.b2j00 = reader.doubles(static_cast<std::size_t>(surfaces));
    validate(result);
    return file;
}

}  // namespace magnetic_coordinate
