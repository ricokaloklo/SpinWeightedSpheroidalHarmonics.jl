# SpinWeightedSpheroidalHarmonics.jl

![license](https://img.shields.io/github/license/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl.svg)](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/SpinWeightedSpheroidalHarmonics.jl)

SpinWeightedSpheroidalHarmonics.jl computes spin-weighted spheroidal harmonics and eigenvalues using a spectral decomposition method. As a natural by-product, it also computes spherical-spheroidal mixing coefficients.

The two main features are implemented as
- `spin_weighted_spheroidal_harmonic` for computing the harmonic ${}_{s} S_{\ell m}(\theta, \phi; c \equiv a \omega)$, and
- `spin_weighted_spheroidal_eigenvalue` for computing the eigenvalue
and both supporting complex spheroidicity $c$ (and hence complex frequency $\omega$). See [Quick-start](@ref) below for some simple examples.

In particular, we use the following normalization convention:
```math
    \int_0^{\pi} \left[ _{s} S_{\ell m}(\theta; c) \right]^{2} \sin(\theta) \, d\theta = \frac{1}{2\pi} \; ,
```
identical to the convention used in the Mathematica package [SpinWeightedSpheroidalHarmonics](https://bhptoolkit.org/SpinWeightedSpheroidalHarmonics/) from the Black Hole Perturbation Toolkit.

Additionally, we provide two similar functions
- `spin_weighted_spherical_harmonic` ${}_{s} Y_{\ell m}(\theta, \phi)$, and
- `spin_weighted_spherical_eigenvalue`
that return the exact harmonic and eigenvalue respectively.

Exact partial derivatives (with respect to either $\theta$ and/or $\phi$) can be evaluated by specifying the derivative order with `theta_derivative` and `phi_derivative` respectively when calling the functions for a harmonic.

Note that v0.5.0 is a *breaking* release.

## Installation
To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("SpinWeightedSpheroidalHarmonics")
```

## Quick-start
### Computing the spin-weighted spheroidal eigenvalue
For example, to compute the spin-weighted spheroidal eigenvalue $\lambda$ for the mode $s = -2, \ell = 2, m = 2, a = 0.7, \omega = 0.5$, simply do
```
using SpinWeightedSpheroidalHarmonics
s=-2; l=2; m=2; a=0.7; omega=0.5;
spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
```

### Computing the spin-weighted spheroidal harmonic
For example, to compute the spin-weighted spheroidal harmonic for the mode $s = -2, \ell = 2, m = 2, a = 0.7, \omega = 0.5$ at $\theta = \pi/6, \phi = \pi/3$, simply do
```
using SpinWeightedSpheroidalHarmonics
s=-2; l=2; m=2; a=0.7; omega=0.5;
theta=π/6; phi=π/3;
# Construct the SpinWeightedSpheroidalHarmonicFunction
swsh = spin_weighted_spheroidal_harmonic(s, l, m, a*omega)
swsh(theta, phi)
```

### Computing the spherical-spheroidal mixing coefficients
For example, to compute the spherical-spheroidal mixing coefficients for the mode $s = -2, \ell = 2, m = 2, a = 0.7, \omega = 0.5$, simply do
```
using SpinWeightedSpheroidalHarmonics
s=-2; l=2; m=2; a=0.7; omega=0.5;
# Construct the SpinWeightedSpheroidalHarmonicFunction
swsh = spin_weighted_spheroidal_harmonic(s, l, m, a*omega)
# swsh.spherical_harmonics_l lists which spherical harmonic \ell modes were used in the decomposition
# swsh.coeffs lists the mixing coefficients for each of the modes
swsh.spherical_harmonics_l, swsh.coeffs
```

## How to cite
If you have used this code in your research that leads to a publication, please cite the following article:
```
@article{Lo:2023fvv,
    author = "Lo, Rico K. L.",
    title = "{Recipes for computing radiation from a Kerr black hole using Generalized Sasaki-Nakamura formalism: I. Homogeneous solutions}",
    eprint = "2306.16469",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "6",
    year = "2023"
}
```

## License
The package is licensed under the MIT License.
