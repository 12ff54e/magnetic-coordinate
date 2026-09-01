#include "boozer_output_internal.hpp"

#include <filesystem>
#include <numbers>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

#include <netcdf.h>

namespace magnetic_coordinate::detail {
namespace {

void check_netcdf(int status, std::string_view operation) {
    if (status != NC_NOERR) {
        throw std::runtime_error("NetCDF " + std::string(operation) + ": " +
                                 nc_strerror(status));
    }
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

void write_boozer_netcdf(const std::filesystem::path& path,
                         const BoozerResult& result,
                         const std::filesystem::path& source_path) {
    const RealBoozerSpectrum spectrum = make_real_spectrum(result);
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    std::error_code ignored;
    std::filesystem::remove(temporary, ignored);

    int file = -1;
    try {
        check_netcdf(nc_create(temporary.string().c_str(), NC_CLOBBER, &file),
                     "create");
        const int surfaces = result.source_ns - result.grid.first_surface;
        const int modes = static_cast<int>(spectrum.m.size());

        int surface_dim = -1;
        int theta_dim = -1;
        int alpha_dim = -1;
        int mode_dim = -1;
        check_netcdf(
            nc_def_dim(file, "surface", static_cast<std::size_t>(surfaces),
                       &surface_dim),
            "define surface dimension");
        check_netcdf(nc_def_dim(file, "theta_b",
                                static_cast<std::size_t>(result.grid.ntheta),
                                &theta_dim),
                     "define theta_b dimension");
        check_netcdf(
            nc_def_dim(file, "alpha",
                       static_cast<std::size_t>(result.grid.nzeta), &alpha_dim),
            "define alpha dimension");
        check_netcdf(nc_def_dim(file, "mode", static_cast<std::size_t>(modes),
                                &mode_dim),
                     "define mode dimension");

        const int surface_dims[1] = {surface_dim};
        const int theta_dims[1] = {theta_dim};
        const int alpha_dims[1] = {alpha_dim};
        const int mode_dims[1] = {mode_dim};
        const int spectrum_dims[2] = {surface_dim, mode_dim};
        const int grid_dims[3] = {surface_dim, alpha_dim, theta_dim};
        int s_var = -1;
        int iota_var = -1;
        int b2j00_var = -1;
        int theta_var = -1;
        int alpha_var = -1;
        int m_var = -1;
        int n_var = -1;
        int rmncc_var = -1;
        int rmnss_var = -1;
        int zmnsc_var = -1;
        int zmncs_var = -1;
        int numnsc_var = -1;
        int numncs_var = -1;
        int b_var = -1;
        int sqrtg_var = -1;
        check_netcdf(nc_def_var(file, "s", NC_DOUBLE, 1, surface_dims, &s_var),
                     "define s");
        check_netcdf(
            nc_def_var(file, "iota", NC_DOUBLE, 1, surface_dims, &iota_var),
            "define iota");
        check_netcdf(
            nc_def_var(file, "b2j00", NC_DOUBLE, 1, surface_dims, &b2j00_var),
            "define b2j00");
        check_netcdf(
            nc_def_var(file, "theta_b", NC_DOUBLE, 1, theta_dims, &theta_var),
            "define theta_b");
        check_netcdf(
            nc_def_var(file, "alpha", NC_DOUBLE, 1, alpha_dims, &alpha_var),
            "define alpha");
        check_netcdf(nc_def_var(file, "mode_m", NC_INT, 1, mode_dims, &m_var),
                     "define mode_m");
        check_netcdf(nc_def_var(file, "mode_n", NC_INT, 1, mode_dims, &n_var),
                     "define mode_n");
        check_netcdf(
            nc_def_var(file, "rmncc", NC_DOUBLE, 2, spectrum_dims, &rmncc_var),
            "define rmncc");
        check_netcdf(
            nc_def_var(file, "rmnss", NC_DOUBLE, 2, spectrum_dims, &rmnss_var),
            "define rmnss");
        check_netcdf(
            nc_def_var(file, "zmnsc", NC_DOUBLE, 2, spectrum_dims, &zmnsc_var),
            "define zmnsc");
        check_netcdf(
            nc_def_var(file, "zmncs", NC_DOUBLE, 2, spectrum_dims, &zmncs_var),
            "define zmncs");
        check_netcdf(nc_def_var(file, "numnsc", NC_DOUBLE, 2, spectrum_dims,
                                &numnsc_var),
                     "define numnsc");
        check_netcdf(nc_def_var(file, "numncs", NC_DOUBLE, 2, spectrum_dims,
                                &numncs_var),
                     "define numncs");
        check_netcdf(nc_def_var(file, "B", NC_DOUBLE, 3, grid_dims, &b_var),
                     "define B");
        check_netcdf(
            nc_def_var(file, "sqrt_g_b", NC_DOUBLE, 3, grid_dims, &sqrtg_var),
            "define sqrt_g_b");

        const auto put_text = [&](const char* name, std::string_view value) {
            check_netcdf(nc_put_att_text(file, NC_GLOBAL, name, value.size(),
                                         value.data()),
                         std::string("write attribute ") + name);
        };
        const auto put_int = [&](const char* name, int value) {
            check_netcdf(
                nc_put_att_int(file, NC_GLOBAL, name, NC_INT, 1, &value),
                std::string("write attribute ") + name);
        };
        const auto put_double = [&](const char* name, double value) {
            check_netcdf(
                nc_put_att_double(file, NC_GLOBAL, name, NC_DOUBLE, 1, &value),
                std::string("write attribute ") + name);
        };
        put_text("schema", BOOZER_SCHEMA);
        put_text("coordinate_convention", COORDINATE_CONVENTION);
        put_text("fourier_convention", FOURIER_CONVENTION);
        put_text("source_path", source_path.string());
        put_int("source_format_version", result.source_format_version);
        put_int("source_ns", result.source_ns);
        put_int("source_ntheta", result.source_ntheta);
        put_int("source_nzeta", result.source_nzeta);
        put_int("source_mpol", result.source_mpol);
        put_int("source_ntor", result.source_ntor);
        put_int("nfp", result.nfp);
        put_int("first_surface", result.grid.first_surface);
        put_int("ntheta", result.grid.ntheta);
        put_int("nzeta", result.grid.nzeta);
        put_int("mmax", result.spectrum.mmax);
        put_int("nmax", result.spectrum.nmax);
        put_int("radial_order", static_cast<int>(result.radial_order));
        put_double("resonance_tolerance", result.resonance_tolerance);

        check_netcdf(nc_enddef(file), "end definitions");
        const auto theta_b = periodic_angles(result.grid.ntheta);
        const auto alpha = periodic_angles(result.grid.nzeta);
        check_netcdf(nc_put_var_double(file, s_var, result.s.data()),
                     "write s");
        check_netcdf(nc_put_var_double(file, iota_var, result.iota.data()),
                     "write iota");
        check_netcdf(
            nc_put_var_double(file, b2j00_var, result.grid.b2j00.data()),
            "write b2j00");
        check_netcdf(nc_put_var_double(file, theta_var, theta_b.data()),
                     "write theta_b");
        check_netcdf(nc_put_var_double(file, alpha_var, alpha.data()),
                     "write alpha");
        check_netcdf(nc_put_var_int(file, m_var, spectrum.m.data()),
                     "write mode_m");
        check_netcdf(nc_put_var_int(file, n_var, spectrum.n.data()),
                     "write mode_n");
        check_netcdf(nc_put_var_double(file, rmncc_var, spectrum.rmncc.data()),
                     "write rmncc");
        check_netcdf(nc_put_var_double(file, rmnss_var, spectrum.rmnss.data()),
                     "write rmnss");
        check_netcdf(nc_put_var_double(file, zmnsc_var, spectrum.zmnsc.data()),
                     "write zmnsc");
        check_netcdf(nc_put_var_double(file, zmncs_var, spectrum.zmncs.data()),
                     "write zmncs");
        check_netcdf(
            nc_put_var_double(file, numnsc_var, spectrum.numnsc.data()),
            "write numnsc");
        check_netcdf(
            nc_put_var_double(file, numncs_var, spectrum.numncs.data()),
            "write numncs");
        check_netcdf(nc_put_var_double(file, b_var, result.grid.b.data()),
                     "write B");
        check_netcdf(
            nc_put_var_double(file, sqrtg_var, result.grid.sqrtg_b.data()),
            "write sqrt_g_b");
        check_netcdf(nc_close(file), "close");
        file = -1;
        std::filesystem::rename(temporary, path);
    } catch (...) {
        if (file >= 0) nc_close(file);
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

}  // namespace magnetic_coordinate::detail
