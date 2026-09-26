# APIs

The main functions are
- [`spin_weighted_spheroidal_eigenvalue`](@ref)
- [`spin_weighted_spheroidal_harmonic`](@ref)
- [`spin_weighted_spherical_eigenvalue`](@ref)
- [`spin_weighted_spherical_harmonic`](@ref)

and their return types use
- [SpectralDecompositionInputParams](@ref)
- [SpinWeightedSpheroidalHarmonicFunction](@ref)
- [SpinWeightedSphericalHarmonicFunction](@ref)

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

## Functions
```@docs
spin_weighted_spheroidal_eigenvalue
```

```@docs
spin_weighted_spheroidal_harmonic
```

```@docs
SpinWeightedSpheroidalHarmonics.SpinWeightedSpheroidalHarmonicFunction
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

## Types
#### SpectralDecompositionInputParams
This is a composite struct type that stores the input parameters for the spectral decomposition 
of a spin-weighted spheroidal harmonic

| field |   |
| :--- | :--- |
| `s` | spin weight $s$ |
| `l` | harmonic index $\ell$, an integer with $\ell \geq \max(|s|, |m|)$ |
| `m` | azimuthal index $m$ |
| `c` | spheroidicity ($c = a\omega$ in the context of BHPT) |
| `N` | number of terms to use in the spectral decomposition (or in the power series, for the `leaver` method) |

#### SpinWeightedSpheroidalHarmonicFunction
This is a composite struct type that stores the output from [`spin_weighted_spheroidal_harmonic`](@ref)

!!! tip

    `SpinWeightedSpheroidalHarmonicFunction(theta, phi)` will return the value of the harmonic at the 
    coordinate $(\theta, \phi)$. For more details, see [`SpinWeightedSpheroidalHarmonics.SpinWeightedSpheroidalHarmonicFunction`](@ref)

| field |   |
| :--- | :--- |
| `params` | a [SpectralDecompositionInputParams](@ref) object storing the input parameters for the spectral decomposition |
| `coeffs` | spectral decomposition coefficients (or the coefficients $a_n$ of Leaver's power series, for the `leaver` method) |
| `spherical_harmonics_l` | an array of [SpinWeightedSphericalHarmonicFunction](@ref) used in the spectral decomposition (empty for the `leaver` method) |
| `normalization_const` | normalization constant to be *divided* to ensure the normalization convention is satisfied |
| `method` | the method used to solve for the harmonic |
| `lambda` | spin-weighted spheroidal eigenvalue $\lambda$ |
| `chebyshev_solution` | numerical solution expressed in Chebyshev polynomials, when solved with the `chebyshev` method |
| `leaver_solution` | Leaver's power-series solution, when solved with the `leaver` method |

#### SpinWeightedSphericalHarmonicFunction
This is a composite struct type that stores information about a spin-weighted spherical harmonic

!!! tip

    `SpinWeightedSphericalHarmonicFunction(theta, phi)` will return the value of the harmonic at the 
    coordinate $(\theta, \phi)$. For more details, see [`SpinWeightedSpheroidalHarmonics.SpinWeightedSphericalHarmonicFunction`](@ref)

| field |   |
| :--- | :--- |
| `s` | spin weight $s$ |
| `l` | harmonic index $\ell$ |
| `m` | azimuthal index $m$ |
| `lambda` | spin-weighted spherical eigenvalue $\lambda$ |
| `method` | the method used to solve for the harmonic |
| `chebyshev_solution` | numerical solution expressed in Chebyshev polynomials |
