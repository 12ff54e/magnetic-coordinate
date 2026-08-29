#include "magnetic_coordinate/boozer_binary.hpp"

#include "magnetic_coordinate/boozer_output.hpp"

#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
#include <netcdf.h>
#endif
#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
#include <hdf5.h>
#endif

#include <cmath>
#include <complex>
#include <exception>
#include <filesystem>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void expect(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

magnetic_coordinate::BoozerResult make_result() {
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
    original.grid.sqrtg_b = {-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0};
    original.grid.b2j00 = {-10.0, -11.0};
    original.spectrum.source_ns = 3;
    original.spectrum.first_surface = 1;
    original.spectrum.mmax = 0;
    original.spectrum.nmax = 0;
    original.spectrum.m = {0};
    original.spectrum.n = {0};
    original.spectrum.r = {{1.0, 0.0}, {3.0, 0.0}};
    original.spectrum.z = {{0.0, 0.0}, {0.0, 0.0}};
    original.spectrum.nu = {{0.0, 0.0}, {0.0, 0.0}};
    return original;
}

void test_round_trip() {
    const auto original = make_result();

    const auto path = std::filesystem::temp_directory_path() /
                      "magnetic-coordinate-boozer-roundtrip.bin";
    magnetic_coordinate::write_boozer_binary(path, original,
                                             "source-output.bin");
    const auto loaded = magnetic_coordinate::read_boozer_binary(path);
    std::filesystem::remove(path);
    expect(loaded.source_path == "source-output.bin", "source provenance");
    expect(loaded.result.source_format_version == 8, "source version");
    expect(loaded.result.nfp == 5, "nfp");
    expect(loaded.result.resonance_tolerance == 2.0e-11, "resonance tolerance");
    expect(loaded.result.s == original.s, "surface coordinates");
    expect(loaded.result.iota == original.iota, "iota");
    expect(loaded.result.spectrum.r == original.spectrum.r, "R spectrum");
    expect(loaded.result.spectrum.z == original.spectrum.z, "Z spectrum");
    expect(loaded.result.spectrum.nu == original.spectrum.nu, "nu spectrum");
    expect(loaded.result.grid.b == original.grid.b, "B grid");
    expect(loaded.result.grid.sqrtg_b == original.grid.sqrtg_b,
           "Jacobian grid");
    expect(loaded.result.grid.b2j00 == original.grid.b2j00, "b2j00");
}

void test_nontrivial_real_parity_round_trip() {
    auto original = make_result();
    original.spectrum.mmax = 1;
    original.spectrum.nmax = 1;
    original.spectrum.m.clear();
    original.spectrum.n.clear();
    for (int n = -1; n <= 1; ++n) {
        for (int m = -1; m <= 1; ++m) {
            original.spectrum.m.push_back(m);
            original.spectrum.n.push_back(n);
        }
    }
    original.spectrum.r.assign(18, {});
    original.spectrum.z.assign(18, {});
    original.spectrum.nu.assign(18, {});
    const auto mode = [](int m, int n) {
        return static_cast<std::size_t>((n + 1) * 3 + (m + 1));
    };
    for (int surface = 0; surface < 2; ++surface) {
        const double scale = std::numbers::pi * std::numbers::pi *
                             static_cast<double>(surface + 1);
        const auto coefficient = [surface](std::size_t mode_index) {
            return static_cast<std::size_t>(surface * 9) + mode_index;
        };
        original.spectrum.r[coefficient(mode(1, 1))] = {1.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(1, -1))] = {2.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(-1, -1))] = {1.5 * scale, 0.0};
        original.spectrum.r[coefficient(mode(-1, 1))] = {2.5 * scale, 0.0};
        original.spectrum.z[coefficient(mode(1, 1))] = {0.0, -0.5 * scale};
        original.spectrum.z[coefficient(mode(1, -1))] = {0.0, -2.0 * scale};
        original.spectrum.z[coefficient(mode(-1, -1))] = {0.0, 0.5 * scale};
        original.spectrum.z[coefficient(mode(-1, 1))] = {0.0, 2.0 * scale};
        original.spectrum.nu[coefficient(mode(1, 1))] = {0.0, -0.5 * scale};
        original.spectrum.nu[coefficient(mode(1, -1))] = {0.0, 1.3 * scale};
        original.spectrum.nu[coefficient(mode(-1, -1))] = {0.0, 0.5 * scale};
        original.spectrum.nu[coefficient(mode(-1, 1))] = {0.0, -1.3 * scale};
    }

    const auto path = std::filesystem::temp_directory_path() /
                      "magnetic-coordinate-boozer-parity.bin";
    magnetic_coordinate::write_boozer_binary(path, original, "source.bin");
    const auto loaded = magnetic_coordinate::read_boozer_binary(path);
    std::filesystem::remove(path);
    const auto close = [](const auto& left, const auto& right) {
        if (left.size() != right.size()) return false;
        for (std::size_t index = 0; index < left.size(); ++index) {
            if (std::abs(left[index] - right[index]) > 1.0e-13) return false;
        }
        return true;
    };
    expect(close(loaded.result.spectrum.r, original.spectrum.r),
           "nontrivial R parity round trip");
    expect(close(loaded.result.spectrum.z, original.spectrum.z),
           "nontrivial Z parity round trip");
    expect(close(loaded.result.spectrum.nu, original.spectrum.nu),
           "nontrivial nu parity round trip");
}

void test_dispatch() {
    using magnetic_coordinate::BoozerOutputFormat;
    expect(
        magnetic_coordinate::resolve_boozer_output_spec("result.BIN").format ==
            BoozerOutputFormat::BINARY,
        "binary suffix dispatch");
    expect(
        magnetic_coordinate::resolve_boozer_output_spec("result.NC").format ==
            BoozerOutputFormat::NETCDF,
        "NetCDF suffix dispatch");
    expect(
        magnetic_coordinate::resolve_boozer_output_spec("result.HDF5").format ==
            BoozerOutputFormat::HDF5,
        "HDF5 suffix dispatch");
#ifndef MAGNETIC_COORDINATE_HAVE_NETCDF
    expect(!magnetic_coordinate::boozer_output_format_available(
               BoozerOutputFormat::NETCDF),
           "NetCDF unavailable without backend");
#endif
#ifndef MAGNETIC_COORDINATE_HAVE_HDF5
    expect(!magnetic_coordinate::boozer_output_format_available(
               BoozerOutputFormat::HDF5),
           "HDF5 unavailable without backend");
#endif
    try {
        static_cast<void>(
            magnetic_coordinate::resolve_boozer_output_spec("result.dat"));
    } catch (const std::invalid_argument&) { return; }
    throw std::runtime_error("unknown output suffix was accepted");
}

#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
void test_netcdf_output() {
    const auto path = std::filesystem::temp_directory_path() /
                      "magnetic-coordinate-boozer-output.nc";
    magnetic_coordinate::write_boozer_output(path, make_result(),
                                             "source-output.bin");
    int file = -1;
    expect(nc_open(path.string().c_str(), NC_NOWRITE, &file) == NC_NOERR,
           "open NetCDF output");
    int variable = -1;
    expect(nc_inq_varid(file, "rmncc", &variable) == NC_NOERR,
           "NetCDF rmncc exists");
    nc_type type = NC_NAT;
    int dimensions = 0;
    expect(nc_inq_var(file, variable, nullptr, &type, &dimensions, nullptr,
                      nullptr) == NC_NOERR &&
               type == NC_DOUBLE && dimensions == 2,
           "NetCDF rmncc is a real rank-2 array");
    expect(nc_inq_varid(file, "R_spectrum", &variable) == NC_ENOTVAR,
           "NetCDF has no complex spectrum variable");
    int complex_dimension = -1;
    expect(nc_inq_dimid(file, "complex_component", &complex_dimension) ==
               NC_EBADDIM,
           "NetCDF has no complex dimension");
    expect(nc_inq_varid(file, "rmncc", &variable) == NC_NOERR,
           "find NetCDF rmncc");
    std::vector<double> rmncc(2);
    expect(nc_get_var_double(file, variable, rmncc.data()) == NC_NOERR,
           "read NetCDF rmncc");
    expect(std::abs(rmncc[0] - 1.0 / (4.0 * std::numbers::pi *
                                      std::numbers::pi)) < 1.0e-15,
           "NetCDF real coefficient normalization");
    expect(nc_close(file) == NC_NOERR, "close NetCDF output");
    std::filesystem::remove(path);
}
#endif

#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
void test_hdf5_output() {
    const auto path = std::filesystem::temp_directory_path() /
                      "magnetic-coordinate-boozer-output.h5";
    magnetic_coordinate::write_boozer_output(path, make_result(),
                                             "source-output.bin");
    const hid_t file =
        H5Fopen(path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    expect(file >= 0, "open HDF5 output");
    const hid_t dataset = H5Dopen2(file, "rmncc", H5P_DEFAULT);
    expect(dataset >= 0, "HDF5 rmncc exists");
    const hid_t space = H5Dget_space(dataset);
    expect(space >= 0 && H5Sget_simple_extent_ndims(space) == 2,
           "HDF5 rmncc is rank 2");
    const hid_t type = H5Dget_type(dataset);
    expect(type >= 0 && H5Tget_class(type) == H5T_FLOAT,
           "HDF5 rmncc is real-valued");
    expect(H5Lexists(file, "R_spectrum", H5P_DEFAULT) == 0,
           "HDF5 has no complex spectrum dataset");
    double rmncc[2] = {};
    expect(H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   rmncc) >= 0,
           "read HDF5 rmncc");
    expect(std::abs(rmncc[0] - 1.0 / (4.0 * std::numbers::pi *
                                      std::numbers::pi)) < 1.0e-15,
           "HDF5 real coefficient normalization");
    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dataset);
    H5Fclose(file);
    std::filesystem::remove(path);
}
#endif

}  // namespace

int main() {
    try {
        test_round_trip();
        test_nontrivial_real_parity_round_trip();
        test_dispatch();
#ifdef MAGNETIC_COORDINATE_HAVE_NETCDF
        test_netcdf_output();
#endif
#ifdef MAGNETIC_COORDINATE_HAVE_HDF5
        test_hdf5_output();
#endif
    } catch (const std::exception& error) {
        std::cerr << "boozer_binary_test: " << error.what() << '\n';
        return 1;
    }
    std::cout << "boozer_binary_test: PASS\n";
    return 0;
}
