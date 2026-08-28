# magnetic-coordinate contributor instructions

- Use strict C++20 and CUDA C++20 without GNU extensions.
- Treat the formulas and normalization in `docs/boozer-binary-v1.md` as
  contracts; change them only with a format-version bump and tests.
- Keep the magnetic axis out of angular root solves.
- Preserve signed Jacobian orientation and apply `nfp` exactly once to
  physical toroidal derivatives.
- GPU FFT operations must be batched across surfaces. CPU B-spline/root work
  may remain surface-local.
- New public inputs require dimension, finiteness, monotonic-map, zero-mode,
  and resonance validation as applicable.
- Add manufactured analytic tests before equilibrium fixtures. Builds use
  warnings as errors; do not suppress project warnings globally.
- Do not commit generated equilibrium, benchmark, plot, or build data.
- Do not put machine names or site-specific absolute paths in tracked files.
