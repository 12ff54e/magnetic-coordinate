cmake_minimum_required(VERSION 3.20)

project(magnetic-coordinate VERSION 0.1.0 LANGUAGES CXX)

include(CTest)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MAGNETIC_COORDINATE_BUILD_TESTING "Build magnetic-coordinate tests" ON)

add_library(magnetic_coordinate STATIC
  src/cumes_binary.cpp
  src/radial.cpp
  src/spectral.cpp)
add_library(magnetic-coordinate::magnetic-coordinate ALIAS magnetic_coordinate)

target_include_directories(magnetic_coordinate
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_compile_features(magnetic_coordinate PUBLIC cxx_std_20)

if(MSVC)
  target_compile_options(magnetic_coordinate PRIVATE /W4 /WX)
else()
  target_compile_options(magnetic_coordinate PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Werror)
endif()

if(BUILD_TESTING AND MAGNETIC_COORDINATE_BUILD_TESTING)
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
endif()
