# SpinWeightedSpheroidalHarmonics.jl

![license](https://img.shields.io/github/license/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl)
[![GitHub release](https://img.shields.io/github/v/release/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl.svg)](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl/releases)
[![Documentation](https://img.shields.io/badge/Documentation-ready)](http://ricokaloklo.github.io/SpinWeightedSpheroidalHarmonics.jl)

SpinWeightedSpheroidalHarmonics.jl computes spin-weighted spheroidal harmonics and eigenvalues using a spectral decomposition method.

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
