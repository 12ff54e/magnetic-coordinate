# Boozer binary format version 1

The file is little-endian. It begins with the eight bytes `MCBOOZ01`, followed
by a signed 32-bit version equal to 1. Strings are a signed 32-bit byte count
followed by that many UTF-8 bytes. All real values are IEEE-754 binary64.
Complex values are stored as consecutive real and imaginary binary64 values.

The header then contains, in order:

1. coordinate-convention string;
2. Fourier-convention string;
3. source output path;
4. source format version, `ns`, `ntheta`, `nzeta`, `mpol`, and `ntor`;
5. `nfp`, `first_surface`, output `ntheta`, and output `nzeta`;
6. spectral `mmax`, `nmax`, radial interpolation order, and resonance
   tolerance.

The coordinate string normatively identifies the grid as
`(s, theta_b, zeta)`: `theta_b` is uniform, while `zeta` is the unchanged
source toroidal coordinate. The full Boozer toroidal coordinate is

```text
zeta_b = zeta + nu.
```

The axis is not exported. Version 1 uses `first_surface=1`, so there are
`ns-first_surface` surface records and
`s = source_surface/(source_ns-1)`.

## Payload

The payload is:

1. `s[surface]`;
2. `iota[surface]`;
3. signed `m[mode]`, then signed `n[mode]`;
4. complex spectra `Rmn[surface][mode]`, `Zmn`, then `numn`;
5. real `B[surface][zeta][theta_b]`;
6. real `sqrtg_b[surface][zeta][theta_b]`;
7. `b2j00[surface]`.

The last index is contiguous. Modes are n-major, with `n=-nmax..nmax` and
`m=-mmax..mmax`. Spectra are continuous Fourier-integral coefficients:

```text
fmn = integral f exp(-i (m theta_b + n zeta)) dtheta_b dzeta.
```

The toroidal integer is a field-period mode. A derivative with respect to the
physical toroidal angle therefore includes `nfp` exactly once. FFTs are
unnormalized internally; the stored coefficient is the raw forward result
times `4*pi^2/(ntheta*nzeta)`.

The stored zero mode obeys

```text
b2j00 = integral B^2 sqrtg_p dtheta_p dzeta,
sqrtg_b = b2j00 / (4*pi^2*B^2).
```

Its sign preserves the source Jacobian orientation.
