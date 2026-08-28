#include "magnetic_coordinate/cumes_binary.hpp"

#include <bit>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename UInt>
void write_little(std::ofstream& output, UInt value) {
    static_assert(std::is_unsigned_v<UInt>);
    for (std::size_t byte = 0; byte < sizeof(UInt); ++byte) {
        const auto value_byte =
            static_cast<char>((value >> (8U * byte)) & static_cast<UInt>(0xff));
        output.write(&value_byte, 1);
    }
}

void write_i32(std::ofstream& output, std::int32_t value) {
    write_little(output, std::bit_cast<std::uint32_t>(value));
}

void write_f64(std::ofstream& output, double value) {
    write_little(output, std::bit_cast<std::uint64_t>(value));
}

void write_u8(std::ofstream& output, std::uint8_t value) {
    output.write(reinterpret_cast<const char*>(&value), 1);
}

void write_string(std::ofstream& output, const std::string& value) {
    write_i32(output, static_cast<std::int32_t>(value.size()));
    output.write(value.data(), static_cast<std::streamsize>(value.size()));
}

void write_f64_vector(std::ofstream& output,
                      const std::vector<double>& values) {
    write_i32(output, static_cast<std::int32_t>(values.size()));
    for (double value : values) write_f64(output, value);
}

void write_empty_vector(std::ofstream& output) {
    write_i32(output, 0);
}

void write_fixture(const std::filesystem::path& path, bool truncate) {
    constexpr int ns = 5;
    constexpr int mpol = 2;
    constexpr int ntor = 1;
    constexpr int mnmax = mpol * (ntor + 1);
    constexpr int ntheta = 3;
    constexpr int nzeta = 2;
    constexpr std::size_t family_size = static_cast<std::size_t>(ns * mnmax);
    constexpr std::size_t half_size =
        static_cast<std::size_t>((ns - 1) * ntheta * nzeta);
    constexpr std::size_t full_size =
        static_cast<std::size_t>(ns * ntheta * nzeta);

    std::ofstream output(path, std::ios::binary);
    output.write("CUMES001", 8);
    write_i32(output, 8);
    write_i32(output, ns);
    write_i32(output, mnmax);
    for (int family = 0; family < 6; ++family) {
        for (std::size_t i = 0; i < family_size; ++i) {
            write_f64(output, static_cast<double>(family * 100) +
                                  static_cast<double>(i));
        }
    }
    write_i32(output, ntheta);
    write_i32(output, nzeta);
    for (int field = 0; field < 7; ++field) {
        for (std::size_t i = 0; i < half_size; ++i) {
            write_f64(output, static_cast<double>(field * 1000) +
                                  static_cast<double>(i));
        }
    }
    for (int field = 0; field < 6; ++field) {
        for (std::size_t i = 0; i < full_size; ++i) {
            write_f64(output, static_cast<double>(field * 2000) +
                                  static_cast<double>(i));
        }
    }
    if (truncate) return;

    write_i32(output, 0);
    write_i32(output, 0);
    write_i32(output, 0);
    write_i32(output, 1);
    write_string(output, "revision");
    write_u8(output, 0);
    for (int string = 0; string < 9; ++string) write_string(output, "");
    write_i32(output, ns);
    write_i32(output, 17);
    write_u8(output, 1);
    write_f64(output, 1.0e-12);
    write_f64(output, 2.0e-12);
    write_f64(output, 3.0e-12);
    write_i32(output, 1);
    write_i32(output, 9);

    write_i32(output, mpol);
    write_i32(output, ntor);
    write_i32(output, 5);
    write_i32(output, ntheta);
    write_i32(output, nzeta);
    write_i32(output, 1);
    write_f64(output, 0.9);
    write_f64(output, -2.5);
    for (int scalar = 0; scalar < 6; ++scalar) write_f64(output, 0.0);
    write_string(output, "cumes-config-v1");
    write_string(output, "power_series");
    write_string(output, "power_series");
    write_string(output, "power_series");
    write_empty_vector(output);
    write_empty_vector(output);
    write_empty_vector(output);
    write_f64_vector(output, {1.0, 0.25});
    write_empty_vector(output);
    write_empty_vector(output);
    write_i32(output, 1);
    write_i32(output, ns);
    write_i32(output, 100);
    write_f64(output, 1.0e-10);
    for (int vector = 0; vector < 10; ++vector) write_empty_vector(output);
    write_i32(output, 0);
    write_i32(output, 1);
    write_string(output, "");
    write_empty_vector(output);
    write_string(output, "");
    write_string(output, "");
    write_i32(output, 0);
}

void test_reads_version_eight() {
    const std::filesystem::path path = "cumes-binary-test.bin";
    write_fixture(path, false);
    const auto equilibrium = magnetic_coordinate::read_cumes_binary(path);
    std::filesystem::remove(path);

    check(equilibrium.format_version == 8, "format version");
    check(equilibrium.ns == 5 && equilibrium.mnmax == 4, "spectral dimensions");
    check(
        equilibrium.mpol == 2 && equilibrium.ntor == 1 && equilibrium.nfp == 5,
        "mode metadata");
    check(equilibrium.ntheta == 3 && equilibrium.nzeta == 2,
          "angular dimensions");
    check(equilibrium.ncurr == 1, "current model");
    check(equilibrium.phiedge == -2.5, "edge flux");
    check(equilibrium.aphi == std::vector<double>({1.0, 0.25}),
          "toroidal flux coefficients");
    check(
        equilibrium.families[magnetic_coordinate::CumesEquilibrium::LMNSC][3] ==
            203.0,
        "lambda spectral family");
    check(equilibrium.half_fields[magnetic_coordinate::CumesEquilibrium::BSUPV]
                                 [4] == 3004.0,
          "half-grid field ordering");
    check(equilibrium.full_fields[0].size() == 30, "full-grid field extent");
    check(equilibrium.native_field_view().sqrtg.size() == 24,
          "native field view");
}

void test_rejects_truncation() {
    const std::filesystem::path path = "cumes-binary-truncated.bin";
    write_fixture(path, true);
    bool rejected = false;
    try {
        static_cast<void>(magnetic_coordinate::read_cumes_binary(path));
    } catch (const std::runtime_error&) { rejected = true; }
    std::filesystem::remove(path);
    check(rejected, "truncated trailer must be rejected");
}

}  // namespace

int main() {
    try {
        test_reads_version_eight();
        test_rejects_truncation();
    } catch (const std::exception& error) {
        std::cerr << "cumes_binary_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "cumes_binary_test: PASS\n";
    return 0;
}
