# Home

**SpinWeightedSpheroidalHarmonics.jl** computes spin-weighted spheroidal harmonics and eigenvalues using a spectral decomposition method.

The two main features are implemented as
- `spin_weighted_spheroidal_harmonic` for computing the harmonic ${}_{s} S_{\ell m}(\theta, \phi; c \equiv a \omega)$, and
- `spin_weighted_spheroidal_eigenvalue` for computing the eigenvalue
and both supporting complex spheroidicity $c$ (and hence frequency $\omega$). See [Quick-start](@ref) below for some simple examples.

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

Internally, there is a caching mechanism that stores the most recent spectral decomposition coefficients such that
when evaluating a harmonic with different values of $\theta$ or $\phi$, no new decomposition will be performed.

## Installation
To install the package using the Julia package manager, simply type the following in the Julia REPL:
```julia
using Pkg
Pkg.add("SpinWeightedSpheroidalHarmonics")
```

## Quick-start
### Computing the spin-weighted spheroidal eigenvalue
For example, to compute the spin-weighted spheroidal eigenvalue $\lambda$ for the mode $s = -2, \ell = 2, m = 2, a = 0.7, \omega = 0.5$, simply do
```@repl
using SpinWeightedSpheroidalHarmonics
s=-2; l=2; m=2; a=0.7; omega=0.5;
spin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)
```

### Computing the spin-weighted spheroidal harmonic
For example, to compute the spin-weighted spheroidal harmonic for the mode $s = -2, \ell = 2, m = 2, a = 0.7, \omega = 0.5$ at $\theta = \pi/6, \phi = \pi/3$, simply do
```@repl
using SpinWeightedSpheroidalHarmonics
s=-2; l=2; m=2; a=0.7; omega=0.5;
theta=π/6; phi=π/3;
spin_weighted_spheroidal_harmonic(s, l, m, a*omega, theta, phi)
```

## License
The package is licensed under the MIT License.