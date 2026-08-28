# magnetic-coordinate

`magnetic-coordinate` converts converged cuMES equilibria into magnetic
coordinate representations. Its first target is a Boozer representation that
keeps the source toroidal angle and exports the shift `nu`, allowing consumers
to construct their own grids uniform in the Boozer toroidal angle.

The implementation is currently being built from the bottom up. The first
available component validates cuMES scientific fields, computes the invariant
magnetic-field strength, unstaggers half-grid quantities onto integer radial
surfaces, and recovers the converged rotational-transform profile.

## Build and test

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Angular B-spline remapping and batched CUDA FFT stages will be added on top of
this host-side preparation layer.
