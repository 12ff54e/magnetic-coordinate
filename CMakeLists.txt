cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED CMAKE_CUDA_HOST_COMPILER AND
   EXISTS "/usr/bin/g++-12")
  set(CMAKE_CUDA_HOST_COMPILER "/usr/bin/g++-12" CACHE FILEPATH
      "CUDA host compiler")
endif()
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 61 CACHE STRING
      "CUDA architectures for magnetic-coordinate")
endif()

project(magnetic-coordinate VERSION 0.2.0 LANGUAGES CXX CUDA)

include(CTest)
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
find_package(CUDAToolkit REQUIRED)

set(_magnetic_coordinate_netcdf_default ON)
set(_magnetic_coordinate_hdf5_default ON)
if(DEFINED CUMES_USE_NETCDF)
  set(_magnetic_coordinate_netcdf_default ${CUMES_USE_NETCDF})
endif()
if(DEFINED CUMES_USE_HDF5)
  set(_magnetic_coordinate_hdf5_default ${CUMES_USE_HDF5})
endif()
option(MAGNETIC_COORDINATE_USE_NETCDF "Enable Boozer NetCDF output"
       ${_magnetic_coordinate_netcdf_default})
option(MAGNETIC_COORDINATE_USE_HDF5 "Enable Boozer HDF5 output"
       ${_magnetic_coordinate_hdf5_default})

set(MAGNETIC_COORDINATE_HAVE_NETCDF FALSE)
set(MAGNETIC_COORDINATE_HAVE_HDF5 FALSE)
set(MAGNETIC_COORDINATE_NETCDF_USES_PKGCONFIG FALSE)
set(MAGNETIC_COORDINATE_HDF5_USES_PKGCONFIG FALSE)
if(MAGNETIC_COORDINATE_USE_NETCDF)
  find_package(netCDF QUIET)
  if(TARGET netCDF::netcdf)
    set(MAGNETIC_COORDINATE_NETCDF_TARGET netCDF::netcdf)
    set(MAGNETIC_COORDINATE_HAVE_NETCDF TRUE)
  elseif(TARGET netCDF::netCDF)
    set(MAGNETIC_COORDINATE_NETCDF_TARGET netCDF::netCDF)
    set(MAGNETIC_COORDINATE_HAVE_NETCDF TRUE)
  else()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
      pkg_check_modules(MAGNETIC_COORDINATE_NETCDF QUIET IMPORTED_TARGET netcdf)
      if(MAGNETIC_COORDINATE_NETCDF_FOUND)
        set(MAGNETIC_COORDINATE_NETCDF_TARGET
            PkgConfig::MAGNETIC_COORDINATE_NETCDF)
        set(MAGNETIC_COORDINATE_NETCDF_USES_PKGCONFIG TRUE)
        set(MAGNETIC_COORDINATE_HAVE_NETCDF TRUE)
      endif()
    endif()
  endif()
endif()
if(MAGNETIC_COORDINATE_USE_HDF5)
  if(TARGET HDF5::HDF5)
    set(MAGNETIC_COORDINATE_HDF5_TARGET HDF5::HDF5)
    set(MAGNETIC_COORDINATE_HAVE_HDF5 TRUE)
  else()
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
      pkg_check_modules(MAGNETIC_COORDINATE_HDF5 QUIET IMPORTED_TARGET hdf5)
      if(MAGNETIC_COORDINATE_HDF5_FOUND)
        set(MAGNETIC_COORDINATE_HDF5_TARGET
            PkgConfig::MAGNETIC_COORDINATE_HDF5)
        set(MAGNETIC_COORDINATE_HDF5_USES_PKGCONFIG TRUE)
        set(MAGNETIC_COORDINATE_HAVE_HDF5 TRUE)
      endif()
    endif()
  endif()
endif()
if(MAGNETIC_COORDINATE_HAVE_NETCDF)
  message(STATUS "magnetic-coordinate: NetCDF output enabled")
else()
  message(STATUS "magnetic-coordinate: NetCDF unavailable; .nc output disabled")
endif()
if(MAGNETIC_COORDINATE_HAVE_HDF5)
  message(STATUS "magnetic-coordinate: HDF5 output enabled")
else()
  message(STATUS "magnetic-coordinate: HDF5 unavailable; .h5 output disabled")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

option(MAGNETIC_COORDINATE_BUILD_TESTING "Build magnetic-coordinate tests" ON)

add_library(magnetic_coordinate STATIC
  src/angular_resample.cpp
  src/boozer_binary.cpp
  src/boozer_output.cpp
  src/boozer_fft.cu
  src/cumes_binary.cpp
  src/mixed_grid.cpp
  src/mixed_spectrum.cu
  src/pest.cpp
  src/radial.cpp
  src/spectral.cpp
  src/spectral_gpu.cu
  src/transform.cpp)
add_library(magnetic-coordinate::magnetic-coordinate ALIAS magnetic_coordinate)
add_library(magnetic-coordinate::magnetic_coordinate ALIAS magnetic_coordinate)

target_include_directories(magnetic_coordinate
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_include_directories(magnetic_coordinate SYSTEM PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/BSplineInterpolation/include)
target_compile_definitions(magnetic_coordinate PRIVATE
  INTP_CELL_LAYOUT
  INTP_PERIODIC_NO_DUMMY_POINT)

target_compile_features(magnetic_coordinate PUBLIC cxx_std_20)
target_link_libraries(magnetic_coordinate PUBLIC CUDA::cudart CUDA::cufft)
if(MAGNETIC_COORDINATE_HAVE_NETCDF)
  target_sources(magnetic_coordinate PRIVATE src/boozer_netcdf.cpp)
  target_compile_definitions(magnetic_coordinate PUBLIC
    MAGNETIC_COORDINATE_HAVE_NETCDF)
  target_link_libraries(magnetic_coordinate PUBLIC
    ${MAGNETIC_COORDINATE_NETCDF_TARGET})
endif()
if(MAGNETIC_COORDINATE_HAVE_HDF5)
  target_sources(magnetic_coordinate PRIVATE src/boozer_hdf5.cpp)
  target_compile_definitions(magnetic_coordinate PUBLIC
    MAGNETIC_COORDINATE_HAVE_HDF5)
  target_link_libraries(magnetic_coordinate PUBLIC
    ${MAGNETIC_COORDINATE_HDF5_TARGET})
endif()
set_target_properties(magnetic_coordinate PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON)

add_executable(cumes-boozer src/main.cpp)
target_link_libraries(cumes-boozer
  PRIVATE magnetic-coordinate::magnetic-coordinate)

configure_package_config_file(
  cmake/magnetic-coordinate-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/magnetic-coordinate-config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/magnetic-coordinate)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/magnetic-coordinate-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(TARGETS magnetic_coordinate cumes-boozer
  EXPORT magnetic-coordinate-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include/magnetic_coordinate
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT magnetic-coordinate-targets
  FILE magnetic-coordinate-targets.cmake
  NAMESPACE magnetic-coordinate::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/magnetic-coordinate)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/magnetic-coordinate-config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/magnetic-coordinate-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/magnetic-coordinate)

if(MSVC)
  target_compile_options(magnetic_coordinate PRIVATE /W4 /WX)
  target_compile_options(cumes-boozer PRIVATE /W4 /WX)
else()
  target_compile_options(magnetic_coordinate PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wpedantic -Wconversion -Werror>
    $<$<COMPILE_LANGUAGE:CUDA>:--Werror cross-execution-space-call>)
  target_compile_options(cumes-boozer PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Werror)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # GCC 12 reports a false maybe-uninitialized warning while destroying the
    # third-party B-spline library's std::variant-based periodic band matrix.
    # Keep every other project and dependency warning under -Werror.
    set_source_files_properties(
      src/angular_resample.cpp src/mixed_grid.cpp src/pest.cpp PROPERTIES
      COMPILE_OPTIONS "-Wno-error=maybe-uninitialized")
  endif()
endif()

if(BUILD_TESTING AND MAGNETIC_COORDINATE_BUILD_TESTING)
  add_executable(magnetic_coordinate_angular_resample_test
    tests/angular_resample_test.cpp)
  target_link_libraries(magnetic_coordinate_angular_resample_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_angular_resample_test
      PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_angular_resample_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.angular-resample
           COMMAND magnetic_coordinate_angular_resample_test)

  add_executable(magnetic_coordinate_boozer_binary_test
    tests/boozer_binary_test.cpp)
  target_link_libraries(magnetic_coordinate_boozer_binary_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_boozer_binary_test
      PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_boozer_binary_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.boozer-binary
           COMMAND magnetic_coordinate_boozer_binary_test)

  add_executable(magnetic_coordinate_radial_test tests/radial_test.cpp)
  target_link_libraries(magnetic_coordinate_radial_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_radial_test PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_radial_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.radial
           COMMAND magnetic_coordinate_radial_test)

  add_executable(magnetic_coordinate_cumes_binary_test
    tests/cumes_binary_test.cpp)
  target_link_libraries(magnetic_coordinate_cumes_binary_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_cumes_binary_test
      PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_cumes_binary_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.cumes-binary
           COMMAND magnetic_coordinate_cumes_binary_test)

  add_executable(magnetic_coordinate_spectral_test tests/spectral_test.cpp)
  target_link_libraries(magnetic_coordinate_spectral_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_spectral_test PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_spectral_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.spectral
           COMMAND magnetic_coordinate_spectral_test)

  add_executable(magnetic_coordinate_spectral_gpu_test
    tests/spectral_gpu_test.cu)
  target_link_libraries(magnetic_coordinate_spectral_gpu_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  set_target_properties(magnetic_coordinate_spectral_gpu_test PROPERTIES
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF)
  add_test(NAME magnetic-coordinate.spectral-gpu
           COMMAND magnetic_coordinate_spectral_gpu_test)

  add_executable(magnetic_coordinate_pest_test tests/pest_test.cpp)
  target_link_libraries(magnetic_coordinate_pest_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_pest_test PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_pest_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.pest
           COMMAND magnetic_coordinate_pest_test)

  add_executable(magnetic_coordinate_mixed_grid_test
    tests/mixed_grid_test.cpp)
  target_link_libraries(magnetic_coordinate_mixed_grid_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  if(MSVC)
    target_compile_options(magnetic_coordinate_mixed_grid_test PRIVATE /W4 /WX)
  else()
    target_compile_options(magnetic_coordinate_mixed_grid_test PRIVATE
      -Wall -Wextra -Wpedantic -Wconversion -Werror)
  endif()
  add_test(NAME magnetic-coordinate.mixed-grid
           COMMAND magnetic_coordinate_mixed_grid_test)

  add_executable(magnetic_coordinate_boozer_fft_test
    tests/boozer_fft_test.cu)
  target_link_libraries(magnetic_coordinate_boozer_fft_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  set_target_properties(magnetic_coordinate_boozer_fft_test PROPERTIES
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF)
  add_test(NAME magnetic-coordinate.boozer-fft
           COMMAND magnetic_coordinate_boozer_fft_test)

  add_executable(magnetic_coordinate_mixed_spectrum_test
    tests/mixed_spectrum_test.cu)
  target_link_libraries(magnetic_coordinate_mixed_spectrum_test
    PRIVATE magnetic-coordinate::magnetic-coordinate)
  set_target_properties(magnetic_coordinate_mixed_spectrum_test PROPERTIES
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF)
  add_test(NAME magnetic-coordinate.mixed-spectrum
           COMMAND magnetic_coordinate_mixed_spectrum_test)
endif()
