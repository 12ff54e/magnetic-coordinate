#include "boozer_output_internal.hpp"

#include <filesystem>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <hdf5.h>

namespace magnetic_coordinate::detail {
namespace {

void check_hdf5(herr_t status, std::string_view operation) {
    if (status < 0) {
        throw std::runtime_error("HDF5 " + std::string(operation) + " failed");
    }
}

void write_string_attribute(hid_t file,
                            const char* name,
                            const std::string& value) {
    const hid_t type = H5Tcopy(H5T_C_S1);
    if (type < 0) throw std::runtime_error("HDF5 copy string type failed");
    const std::size_t width = value.size() + 1;
    hid_t space = -1;
    hid_t attribute = -1;
    try {
        check_hdf5(H5Tset_size(type, width), "set string width");
        space = H5Screate(H5S_SCALAR);
        if (space < 0)
            throw std::runtime_error("HDF5 create scalar space failed");
        attribute =
            H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
        if (attribute < 0) {
            throw std::runtime_error("HDF5 create string attribute failed");
        }
        check_hdf5(H5Awrite(attribute, type, value.c_str()),
                   "write string attribute");
        H5Aclose(attribute);
        H5Sclose(space);
        H5Tclose(type);
    } catch (...) {
        if (attribute >= 0) H5Aclose(attribute);
        if (space >= 0) H5Sclose(space);
        H5Tclose(type);
        throw;
    }
}

template <typename T>
void write_scalar_attribute(hid_t file,
                            const char* name,
                            hid_t type,
                            const T& value) {
    const hid_t space = H5Screate(H5S_SCALAR);
    if (space < 0) throw std::runtime_error("HDF5 create scalar space failed");
    const hid_t attribute =
        H5Acreate2(file, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    if (attribute < 0) {
        throw std::runtime_error("HDF5 create scalar attribute failed");
    }
    const herr_t status = H5Awrite(attribute, type, &value);
    H5Aclose(attribute);
    check_hdf5(status, "write scalar attribute");
}

void write_dataset(hid_t file,
                   const char* name,
                   hid_t type,
                   std::span<const hsize_t> dimensions,
                   const void* data) {
    const hid_t space = H5Screate_simple(static_cast<int>(dimensions.size()),
                                         dimensions.data(), nullptr);
    if (space < 0) throw std::runtime_error("HDF5 create dataset space failed");
    const hid_t dataset = H5Dcreate2(file, name, type, space, H5P_DEFAULT,
                                     H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space);
    if (dataset < 0) throw std::runtime_error("HDF5 create dataset failed");
    const herr_t status =
        H5Dwrite(dataset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dataset);
    check_hdf5(status, "write dataset");
}

std::vector<double> periodic_angles(int count) {
    std::vector<double> values(static_cast<std::size_t>(count));
    for (int index = 0; index < count; ++index) {
        values[static_cast<std::size_t>(index)] = 2.0 * std::numbers::pi *
                                                  static_cast<double>(index) /
                                                  static_cast<double>(count);
    }
    return values;
}

}  // namespace

void write_boozer_hdf5(const std::filesystem::path& path,
                       const BoozerResult& result,
                       const std::filesystem::path& source_path) {
    const RealBoozerSpectrum spectrum = make_real_spectrum(result);
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);
    hid_t file = -1;
    try {
        file = H5Fcreate(temporary.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
        if (file < 0) throw std::runtime_error("HDF5 create failed");

        write_string_attribute(file, "schema", std::string(BOOZER_SCHEMA));
        write_string_attribute(file, "coordinate_convention",
                               std::string(COORDINATE_CONVENTION));
        write_string_attribute(file, "fourier_convention",
                               std::string(FOURIER_CONVENTION));
        write_string_attribute(file, "source_path", source_path.string());
        write_scalar_attribute(file, "source_format_version", H5T_NATIVE_INT,
                               result.source_format_version);
        write_scalar_attribute(file, "source_ns", H5T_NATIVE_INT,
                               result.source_ns);
        write_scalar_attribute(file, "source_ntheta", H5T_NATIVE_INT,
                               result.source_ntheta);
        write_scalar_attribute(file, "source_nzeta", H5T_NATIVE_INT,
                               result.source_nzeta);
        write_scalar_attribute(file, "source_mpol", H5T_NATIVE_INT,
                               result.source_mpol);
        write_scalar_attribute(file, "source_ntor", H5T_NATIVE_INT,
                               result.source_ntor);
        write_scalar_attribute(file, "nfp", H5T_NATIVE_INT, result.nfp);
        write_scalar_attribute(file, "first_surface", H5T_NATIVE_INT,
                               result.grid.first_surface);
        write_scalar_attribute(file, "ntheta", H5T_NATIVE_INT,
                               result.grid.ntheta);
        write_scalar_attribute(file, "nzeta", H5T_NATIVE_INT,
                               result.grid.nzeta);
        write_scalar_attribute(file, "mmax", H5T_NATIVE_INT,
                               result.spectrum.mmax);
        write_scalar_attribute(file, "nmax", H5T_NATIVE_INT,
                               result.spectrum.nmax);
        const int radial_order = static_cast<int>(result.radial_order);
        write_scalar_attribute(file, "radial_order", H5T_NATIVE_INT,
                               radial_order);
        write_scalar_attribute(file, "resonance_tolerance", H5T_NATIVE_DOUBLE,
                               result.resonance_tolerance);

        const hsize_t surface_dims[1] = {static_cast<hsize_t>(result.s.size())};
        const hsize_t theta_dims[1] = {
            static_cast<hsize_t>(result.grid.ntheta)};
        const hsize_t zeta_dims[1] = {static_cast<hsize_t>(result.grid.nzeta)};
        const hsize_t mode_dims[1] = {static_cast<hsize_t>(spectrum.m.size())};
        const hsize_t spectrum_dims[2] = {
            static_cast<hsize_t>(result.s.size()),
            static_cast<hsize_t>(spectrum.m.size())};
        const hsize_t grid_dims[3] = {static_cast<hsize_t>(result.s.size()),
                                      static_cast<hsize_t>(result.grid.nzeta),
                                      static_cast<hsize_t>(result.grid.ntheta)};
        const auto theta_b = periodic_angles(result.grid.ntheta);
        const auto zeta = periodic_angles(result.grid.nzeta);
        write_dataset(file, "s", H5T_NATIVE_DOUBLE, surface_dims,
                      result.s.data());
        write_dataset(file, "iota", H5T_NATIVE_DOUBLE, surface_dims,
                      result.iota.data());
        write_dataset(file, "b2j00", H5T_NATIVE_DOUBLE, surface_dims,
                      result.grid.b2j00.data());
        write_dataset(file, "theta_b", H5T_NATIVE_DOUBLE, theta_dims,
                      theta_b.data());
        write_dataset(file, "zeta", H5T_NATIVE_DOUBLE, zeta_dims, zeta.data());
        write_dataset(file, "mode_m", H5T_NATIVE_INT, mode_dims,
                      spectrum.m.data());
        write_dataset(file, "mode_n", H5T_NATIVE_INT, mode_dims,
                      spectrum.n.data());
        write_dataset(file, "rmncc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.rmncc.data());
        write_dataset(file, "rmnss", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.rmnss.data());
        write_dataset(file, "zmnsc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.zmnsc.data());
        write_dataset(file, "zmncs", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.zmncs.data());
        write_dataset(file, "numnsc", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.numnsc.data());
        write_dataset(file, "numncs", H5T_NATIVE_DOUBLE, spectrum_dims,
                      spectrum.numncs.data());
        write_dataset(file, "B", H5T_NATIVE_DOUBLE, grid_dims,
                      result.grid.b.data());
        write_dataset(file, "sqrt_g_b", H5T_NATIVE_DOUBLE, grid_dims,
                      result.grid.sqrtg_b.data());
        check_hdf5(H5Fclose(file), "close");
        file = -1;
        std::filesystem::rename(temporary, path);
    } catch (...) {
        if (file >= 0) H5Fclose(file);
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

}  // namespace magnetic_coordinate::detail
