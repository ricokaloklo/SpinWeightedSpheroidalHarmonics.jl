# APIs

The main functions are
- [`spin_weighted_spheroidal_eigenvalue`](@ref)
- [`spin_weighted_spheroidal_harmonic`](@ref)
- [`spin_weighted_spherical_eigenvalue`](@ref)
- [`spin_weighted_spherical_harmonic`](@ref)
- [`value_and_theta_derivative`](@ref), which evaluates a harmonic and its ``\theta``-derivative in one pass

and their return types use
- [SpectralDecompositionInputParams](@ref)
- [SpinWeightedSpheroidalHarmonicFunction](@ref)
- [SpinWeightedSphericalHarmonicFunction](@ref)

Every function computes in the floating-point type of its input, and the functions that take a
spheroidicity also accept the keyword `precision`; see [Arbitrary precision](@ref).

## Complex frequencies

```julia
c = 1.0 - 4.0im
lambda = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, c)
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, c)
S(1.1, 0.3)
S(1.1, 0.3; theta_derivative=1)
```

To follow a mode along an ordered path, pass the complex points to
`track_angular_mode`. The returned `states` include intermediate points.

```julia
path = [0.0, 0.5 - 2.0im, 1.0 - 4.0im]
result = track_angular_mode(-2, 2, 2, path)
pair = last(result.states)
lambda = pair.lambda
S = spin_weighted_spheroidal_harmonic(pair)
```

Use `backend=:dense` with `spin_weighted_spheroidal_harmonic` or
`track_angular_mode` to compare with the default backend along the same path.

For complex `c` the harmonic is normalized with its modulus,
``\int_0^\pi |{}_s S_{\ell m}|^2 \sin\theta \, d\theta = 1/(2\pi)``, and every method fixes its overall
phase by making the overlap with ``{}_s Y_{\ell m}`` real and positive. Because of the complex
conjugate, the normalization constant, and with it the harmonic, is not an analytic function of `c`.

With `normalization=:analytic` the harmonic is instead normalized with
``\int_0^\pi \left[{}_s S_{\ell m}(\theta, 0)\right]^2 \sin\theta \, d\theta = 1/(2\pi)``, without the
complex conjugate, as the Mathematica package of the Black Hole Perturbation Toolkit does. The harmonic
is then the analytic continuation in `c` of the harmonic for real `c`, with its overlap with
``{}_s Y_{\ell m}`` the continuation of a positive number (no longer real). Use it for quasinormal-mode
residues and anything else that continues the harmonic to complex ``c = a\omega``:

```julia
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, 0.5 - 0.5im; normalization=:analytic)
S = spin_weighted_spheroidal_harmonic(pair; normalization=:analytic)  # for a continued eigenpair
```

The two normalizations agree for real `c`, and the eigenvalue does not depend on the normalization.
The default stays `normalization=:modulus`.

## Functions
```@docs
spin_weighted_spheroidal_eigenvalue
```

```@docs
spin_weighted_spheroidal_harmonic(::Int, ::Any, ::Int, ::Any)
```

```@docs
SpinWeightedSpheroidalHarmonics.SpinWeightedSpheroidalHarmonicFunction
```

```@docs
value_and_theta_derivative
```

```@docs
spin_weighted_spherical_eigenvalue
```

```@docs
spin_weighted_spherical_harmonic
```

```@docs
SpinWeightedSpheroidalHarmonics.SpinWeightedSphericalHarmonicFunction
```

## Continuation in complex spheroidicity

```@docs
continue_angular_mode
track_angular_mode
spin_weighted_spheroidal_harmonic(::AngularEigenpair)
angular_observables
continue_angular_lateral_pair
angular_precision_certificate
mirror_angular_parameters
angular_mirror_residual
AngularEigenpair
AngularPathResult
AngularContinuationError
AngularCache
clear_angular_cache!
```

## Types
#### SpectralDecompositionInputParams
This is a composite struct type that stores the input parameters for the spectral decomposition 
of a spin-weighted spheroidal harmonic

| field |   |
| :--- | :--- |
| `s` | spin weight ``s`` |
| `l` | harmonic index ``\ell``, an integer with ``\ell \geq \max(\lvert s \rvert, \lvert m \rvert)`` |
| `m` | azimuthal index ``m`` |
| `c` | spheroidicity (``c = a\omega`` in the context of BHPT) |
| `N` | number of terms to use in the spectral decomposition (or in the power series, for the `leaver` method) |

#### SpinWeightedSpheroidalHarmonicFunction
This is a composite struct type that stores the output from [`spin_weighted_spheroidal_harmonic`](@ref)

!!! tip

    `SpinWeightedSpheroidalHarmonicFunction(theta, phi)` will return the value of the harmonic at the 
    coordinate ``(\theta, \phi)``. For more details, see [`SpinWeightedSpheroidalHarmonics.SpinWeightedSpheroidalHarmonicFunction`](@ref)

| field |   |
| :--- | :--- |
| `params` | a [SpectralDecompositionInputParams](@ref) object storing the input parameters for the spectral decomposition |
| `coeffs` | spectral decomposition coefficients (or the coefficients ``a_n`` of Leaver's power series, for the `leaver` method) |
| `spherical_harmonics_l` | an array of [SpinWeightedSphericalHarmonicFunction](@ref) used in the spectral decomposition (empty for the `leaver` method) |
| `normalization_const` | normalization constant to be *divided* to ensure the normalization convention is satisfied |
| `method` | the method used to solve for the harmonic |
| `lambda` | spin-weighted spheroidal eigenvalue ``\lambda`` |
| `chebyshev_solution` | numerical solution expressed in Chebyshev polynomials, when solved with the `chebyshev` method |
| `leaver_solution` | Leaver's power-series solution, when solved with the `leaver` method |
| `precision` | the `BigFloat` precision the harmonic was built with, at which it is evaluated (`0` for `Float64`); the type parameter is the floating-point type it was built in |

#### SpinWeightedSphericalHarmonicFunction
This is a composite struct type that stores information about a spin-weighted spherical harmonic

!!! tip

    `SpinWeightedSphericalHarmonicFunction(theta, phi)` will return the value of the harmonic at the 
    coordinate ``(\theta, \phi)``. For more details, see [`SpinWeightedSpheroidalHarmonics.SpinWeightedSphericalHarmonicFunction`](@ref)

| field |   |
| :--- | :--- |
| `s` | spin weight ``s`` |
| `l` | harmonic index ``\ell`` |
| `m` | azimuthal index ``m`` |
| `lambda` | spin-weighted spherical eigenvalue ``\lambda`` |
| `method` | the method used to solve for the harmonic |
| `chebyshev_solution` | numerical solution expressed in Chebyshev polynomials |
| `jacobi_theta_terms` | the first and second theta-derivative terms, when evaluated with the `jacobi` method |
| `direct_terms` | the prepared terms of the direct sum, when evaluated with the `direct` method |
| `precision` | the `BigFloat` precision the harmonic was built with, at which it is evaluated (`0` for `Float64`); the type parameter is the floating-point type it was built in |
