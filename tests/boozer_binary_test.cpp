#include "magnetic_coordinate/boozer_binary.hpp"

#include <complex>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void expect(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void test_round_trip() {
    magnetic_coordinate::BoozerResult original;
    original.source_format_version = 8;
    original.source_ns = 3;
    original.source_ntheta = 8;
    original.source_nzeta = 6;
    original.source_mpol = 3;
    original.source_ntor = 2;
    original.nfp = 5;
    original.radial_order =
        magnetic_coordinate::RadialInterpolationOrder::TWO_POINT;
    original.resonance_tolerance = 2.0e-11;
    original.s = {0.5, 1.0};
    original.iota = {0.4, 0.5};
    original.grid.source_ns = 3;
    original.grid.first_surface = 1;
    original.grid.ntheta = 2;
    original.grid.nzeta = 2;
    original.grid.b = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    original.grid.sqrtg_b = {-1.0, -2.0, -3.0, -4.0,
                             -5.0, -6.0, -7.0, -8.0};
    original.grid.b2j00 = {-10.0, -11.0};
    original.spectrum.source_ns = 3;
    original.spectrum.first_surface = 1;
    original.spectrum.mmax = 0;
    original.spectrum.nmax = 0;
    original.spectrum.m = {0};
    original.spectrum.n = {0};
    original.spectrum.r = {{1.0, 2.0}, {3.0, 4.0}};
    original.spectrum.z = {{5.0, 6.0}, {7.0, 8.0}};
    original.spectrum.nu = {{0.0, 0.0}, {0.0, 0.0}};

    const auto path = std::filesystem::temp_directory_path() /
                      "magnetic-coordinate-boozer-roundtrip.bin";
    magnetic_coordinate::write_boozer_binary(path, original,
                                               "source-output.bin");
    const auto loaded = magnetic_coordinate::read_boozer_binary(path);
    std::filesystem::remove(path);
    expect(loaded.source_path == "source-output.bin", "source provenance");
    expect(loaded.result.source_format_version == 8, "source version");
    expect(loaded.result.nfp == 5, "nfp");
    expect(loaded.result.resonance_tolerance == 2.0e-11,
           "resonance tolerance");
    expect(loaded.result.s == original.s, "surface coordinates");
    expect(loaded.result.iota == original.iota, "iota");
    expect(loaded.result.spectrum.r == original.spectrum.r, "R spectrum");
    expect(loaded.result.spectrum.z == original.spectrum.z, "Z spectrum");
    expect(loaded.result.spectrum.nu == original.spectrum.nu,
           "nu spectrum");
    expect(loaded.result.grid.b == original.grid.b, "B grid");
    expect(loaded.result.grid.sqrtg_b == original.grid.sqrtg_b,
           "Jacobian grid");
    expect(loaded.result.grid.b2j00 == original.grid.b2j00, "b2j00");
}

}  // namespace

int main() {
    try {
        test_round_trip();
    } catch (const std::exception& error) {
        std::cerr << "boozer_binary_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "boozer_binary_test: PASS\n";
    return 0;
}
