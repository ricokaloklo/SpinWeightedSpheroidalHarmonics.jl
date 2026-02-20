# APIs

There are 4 functions that are exported, namely
- [`spin_weighted_spheroidal_eigenvalue`](@ref)
- [`spin_weighted_spheroidal_harmonic`](@ref)
- [`spin_weighted_spherical_eigenvalue`](@ref)
- [`spin_weighted_spherical_harmonic`](@ref)

and there are 3 custom types that are exported, i.e.
- [SpectralDecompositionInputParams](@ref)
- [SpinWeightedSpheroidalHarmonicFunction](@ref)
- [SpinWeightedSphericalHarmonicFunction](@ref)

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
| `l` | harmonic index $\ell$ |
| `m` | azimuthal index $m$ |
| `c` | spheroidicity ($c = a\omega$ in the context of BHPT) |
| `N` | number of terms to use in the spectral decomposition |

#### SpinWeightedSpheroidalHarmonicFunction
This is a composite struct type that stores the output from [`spin_weighted_spheroidal_harmonic`](@ref)

!!! tip

    `SpinWeightedSpheroidalHarmonicFunction(theta, phi)` will return the value of the harmonic at the 
    coordinate $(\theta, \phi)$. For more details, see [`SpinWeightedSpheroidalHarmonics.SpinWeightedSpheroidalHarmonicFunction`](@ref)

| field |   |
| :--- | :--- |
| `params` | a [SpectralDecompositionInputParams](@ref) object storing the input parameters for the spectral decomposition |
| `coeffs` | spectral decomposition coefficients |
| `spherical_harmonics_l` | an array of [SpinWeightedSphericalHarmonicFunction](@ref) used in the spectral decomposition |
| `normalization_const` | normalization constant to be *divided* to ensure the normalization convention is satisfied |
| `method` | the method used to solve for the harmonic |
| `lambda` | spin-weighted spheroidal eigenvalue $\lambda$ |

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
