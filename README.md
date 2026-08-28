# magnetic-coordinate

`magnetic-coordinate` converts a converged cuMES schema-v8 equilibrium into a
Boozer representation. It exports geometry on a grid uniform in Boozer
poloidal angle while deliberately retaining the source toroidal angle:

```text
(s, theta_b, zeta),       zeta_b = zeta + nu(theta_b, zeta)
```

This is a mixed grid, not a grid uniform in both Boozer angles. Spectral `nu`
is included so consumers can construct whatever uniform `zeta_b` grid they
need.

The transform validates cuMES scientific fields, computes invariant `B`,
unstaggers radial fields, recovers converged iota, reconstructs physical PEST
lambda, performs both monotone angular inversions with periodic cubic
B-splines, solves the magnetic differential equation with batched CUDA FFTs,
and analyzes the final mixed-grid `R`, `Z`, and `nu` spectra on the GPU.

## Build and test

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure

./build/cumes-boozer cumes-output.bin --output boozer-output.bin
```

The default output poloidal resolution and spectral truncation follow the
source equilibrium. They can be selected independently with `--ntheta`,
`--nzeta`, `--mmax`, and `--nmax`; changing `nzeta` periodically resamples the
same unchanged zeta coordinate. `--radial-order 2|4` controls only half-to-full
radial interpolation. See [the format contract](docs/boozer-binary-v1.md) for
exact array ordering and Fourier normalization.

The magnetic axis is intentionally excluded because flux angles are
degenerate there; the file records `first_surface=1` explicitly.

## Install and consume

```bash
cmake --install build --prefix /desired/prefix
```

The install exports `magnetic-coordinate::magnetic_coordinate` and the
`cumes-boozer` executable. A consuming CMake project can use
`find_package(magnetic-coordinate CONFIG REQUIRED)` and link the exported
target.

The same transform is available without an intermediate file. Populate a
non-owning `magnetic_coordinate::CumesEquilibriumView` with spans into a
solver's converged arrays, then call:

```cpp
#include <magnetic_coordinate/boozer_binary.hpp>
#include <magnetic_coordinate/transform.hpp>

magnetic_coordinate::TransformSettings settings;
auto result = magnetic_coordinate::transform_to_boozer(view, settings);

// Optional: keep `result` in memory, or serialize it explicitly.
magnetic_coordinate::write_boozer_binary("boozer-output.bin", result,
                                         "in-memory equilibrium");
```

`CumesEquilibriumView` does not copy or own its arrays; they must remain valid
until `transform_to_boozer` returns. `transform_cumes_file` and the standalone
executable adapt file-backed storage to this same interface, so both paths run
the identical implementation.
