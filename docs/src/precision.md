# Arbitrary precision

By default every computation is in `Float64`. The package also computes in `BigFloat`, at any
number of bits: the precision is set by the input, and nothing else changes in how a harmonic or an
eigenvalue is computed or used.

## Usage

There are three ways to ask for more digits.

**`BigFloat` input.** Every result is computed in the floating-point type of `c` (and, for a
harmonic, evaluated in the wider of that type and the type of the angles):

```julia
setprecision(BigFloat, 256) do
    spin_weighted_spheroidal_eigenvalue(-2, 2, 2, big"0.35")
end
```

**The `precision` keyword.** `precision = p` converts the input to `p`-bit `BigFloat`:

```@example precision
using SpinWeightedSpheroidalHarmonics

λ64 = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.35)           # Float64, the default
λ = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.35; precision=256)
```

```@example precision
λ - λ64     # the Float64 error, since BigFloat(0.35) is the Float64 value of 0.35
```

**Rational input.** With `precision`, a rational number is rounded once, in the target
precision, so `7//20` is ``7/20`` to the last of the 256 bits:

```@example precision
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, 7//20; precision=256)
S(big(11) / 10, 0)
```

Complex `c` is supported in the same way (`precision` converts the real and imaginary parts).
The keyword is accepted by `spin_weighted_spheroidal_harmonic`, `spin_weighted_spheroidal_eigenvalue`,
`spin_weighted_spherical_harmonic`, `continue_angular_mode`, `track_angular_mode` and
`continue_angular_lateral_pair`. Every other function follows the type of its arguments (a harmonic
built from an `AngularEigenpair{BigFloat}` is a `BigFloat` harmonic). `spin_weighted_spherical_eigenvalue`
returns the exact integer.

Every method is available in `BigFloat`: the spectral decomposition (`"auto"`, `"direct"`, `"jacobi"`,
both backends), `"chebyshev"` and `"leaver"`.

## A harmonic remembers its precision

A harmonic evaluates at the precision it was built with, whatever `setprecision` is in effect where
it is called:

```@example precision
v = setprecision(BigFloat, 64) do
    S(1.1, 0.0; theta_derivative=1)
end
precision(real(v))
```

Harmonics, eigenvalues and continued modes can be computed and evaluated from several threads at
once, at any precision. The only shared state is the cache of continued modes
(`DEFAULT_ANGULAR_CACHE`), which is locked and keyed by everything a result depends on, including
the number type and its precision; `clear_angular_cache!` empties it.

## What follows the precision

The eigenpairs of the spectral decomposition are found in `Float64` with LAPACK and refined in the
working precision on the band matrix (Rayleigh quotient iteration); in `Float64` the LAPACK result
is returned directly. Where two modes are too close for `Float64` to tell them apart, the mode is
selected in the working precision (by its position among the eigenvalues for real `c`, by following
the path for complex `c`), so a wider precision also resolves harmonics that `Float64` cannot. The
matrix size, Leaver's continued fraction and power series, the Chebyshev representations and the
tolerances all follow the precision, which keeps each tolerance at its place between the rounding
and the physical scale. At 256 bits eigenvalues and harmonics are accurate to about ``10^{-70}``,
less only where the mode itself is ill-conditioned (nearly degenerate with a neighbour).

## Accuracy

Relative errors against a 512-bit reference, the largest over ``\theta`` from pole to pole, for the
eigenvalue, and for the harmonic with its first two ``\theta``-derivatives relative to the size of each
(default method):

| mode, spheroidicity | | `Float64` | 128 bits | 256 bits |
| :--- | :--- | :--- | :--- | :--- |
| ``(-2, 2, 2)``, ``a = 0.7``, ``\omega = 0.5326 - 0.0808i`` | ``\lambda`` | ``3 \times 10^{-16}`` | ``3 \times 10^{-39}`` | ``2 \times 10^{-77}`` |
| | ``S, \partial_\theta S, \partial_\theta^2 S`` | ``2 \times 10^{-15}`` | ``2 \times 10^{-38}`` | ``4 \times 10^{-76}`` |
| ``(-2, 3, 1)``, ``c = 1 - 2i`` | ``\lambda`` | ``2 \times 10^{-14}`` | ``1 \times 10^{-38}`` | ``5 \times 10^{-77}`` |
| | ``S, \partial_\theta S, \partial_\theta^2 S`` | ``2 \times 10^{-14}`` | ``1 \times 10^{-37}`` | ``2 \times 10^{-76}`` |
| ``(-2, 20, 2)``, ``c = -1.5 - 0.7i`` | ``\lambda`` | ``2 \times 10^{-16}`` | ``1 \times 10^{-39}`` | ``3 \times 10^{-77}`` |
| | ``S, \partial_\theta S, \partial_\theta^2 S`` | ``9 \times 10^{-14}`` | ``2 \times 10^{-37}`` | ``1 \times 10^{-75}`` |

The derivatives keep this accuracy at the poles, for every method except the Chebyshev solver
(`method="chebyshev"`), which loses a few digits in its derivatives near ``\theta = 0`` and ``\pi``.

## Modes that `Float64` cannot resolve

When a mode is closer to a neighbour than `Float64` can tell apart, the harmonic is not returned
in `Float64`: an error names the mode and the gap. The same mode is resolved in a wider precision,
given in that precision:

```julia
spin_weighted_spheroidal_harmonic(-2, 2, 2, -20.0)                    # error: nearly degenerate
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, -20; precision=256)  # resolved
```

## Three points about the input

- `BigFloat(0.35)` is the binary value of the `Float64` number `0.35`, not ``7/20``. To give a
  decimal value, pass a rational (`7//20`) with `precision`, or `parse(BigFloat, "0.35")`.
- `big"0.35"` is parsed when the code is compiled, at the precision in effect then. Inside a
  function, use `parse(BigFloat, "0.35")` or a rational.
- A near-degeneracy (two modes closer than the rounding) is judged in the precision of the input.
  A spheroidicity that is unresolvable in `Float64` can be resolvable at 256 bits; pass it in the
  target precision.

## Cost

`BigFloat` is roughly a thousand times slower than `Float64`. At 256 bits, after compilation, an
eigenvalue takes about 0.05 s, building a harmonic 0.01–0.5 s (Leaver's method a few seconds for
``|c| \gtrsim 4``, Chebyshev about 0.5 s) and an evaluation about 0.02 s for the spectral
decomposition; the cost grows with the number of bits and with ``|c|`` and ``\ell``.

## Supported types

`Float64` (the default) and `BigFloat`. `Float32`, integer and rational inputs are computed in
`Float64`. Other formats (double-double, quadruple, multi-word floats) are not supported.
