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

project(magnetic-coordinate VERSION 0.1.0 LANGUAGES CXX CUDA)

include(CTest)
find_package(CUDAToolkit REQUIRED)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_EXTENSIONS OFF)

option(MAGNETIC_COORDINATE_BUILD_TESTING "Build magnetic-coordinate tests" ON)

add_library(magnetic_coordinate STATIC
  src/boozer_binary.cpp
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
set_target_properties(magnetic_coordinate PROPERTIES
  CUDA_SEPARABLE_COMPILATION ON)

add_executable(cumes-boozer src/main.cpp)
target_link_libraries(cumes-boozer
  PRIVATE magnetic-coordinate::magnetic-coordinate)

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
    set_source_files_properties(src/mixed_grid.cpp src/pest.cpp PROPERTIES
      COMPILE_OPTIONS "-Wno-error=maybe-uninitialized")
  endif()
endif()

if(BUILD_TESTING AND MAGNETIC_COORDINATE_BUILD_TESTING)
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
