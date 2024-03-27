module SpinWeightedSpheroidalHarmonics

using LinearAlgebra
using QuadGK

include("harmonic.jl")
include("spectral.jl")

export spin_weighted_spheroidal_harmonic, spin_weighted_spherical_harmonic, spin_weighted_spheroidal_eigenvalue, spin_weighted_spherical_eigenvalue # Expose these functions to the user
export Teukolsky_lambda_const # For backward compatbility

struct SpectralDecompositionInputParams
    s::Int
    l::Int
    m::Int
    c
    N::Int
end

struct SpinWeightedSpheroidalHarmonicFunction
    params::SpectralDecompositionInputParams
    coeffs
    normalization_const
    lambda
end

# Implement pretty printing for SpinWeightedSpheroidalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSpheroidalHarmonicFunction)
    print(io, "SpinWeightedSpheroidalHarmonicFunction(s = $(swsh_func.params.s), l = $(swsh_func.params.l), m = $(swsh_func.params.m), c = $(swsh_func.params.c), lambda = $(swsh_func.lambda))")
end

# Not really necessary but this is to maintain uniformity
struct SpinWeightedSphericalHarmonicFunction
    s::Int
    l::Int
    m::Int
    lambda
end

# Implement pretty printing for SpinWeightedSphericalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSphericalHarmonicFunction)
    print(io, "SpinWeightedSphericalHarmonicFunction(s = $(swsh_func.s), l = $(swsh_func.l), m = $(swsh_func.m), lambda = $(swsh_func.lambda))")
end

function _unnormalized_spin_weighted_spheroidal_harmonic(coefficients_params, coefficients, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)
    # Compute the spin-weighted spherical harmonics needed
    output = 0.0
    spherical_harmonics_l = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)

    for (idx, spherical_harmonic_l) in enumerate(spherical_harmonics_l)
        output += coefficients[idx] * _nth_derivative_spherical_harmonic(coefficients_params.s, spherical_harmonic_l, coefficients_params.m, theta_derivative, phi_derivative, theta, phi)
    end
    output
end

function _determine_matrix_size_N(s::Int, l::Int, m::Int)
    #=
    Determine a suitable value of N for the spectral decomposition

    The value of N calculated here is essentially lmax for
    the spectral decomposition. Then we apply a 'buffer' of 10
    =#
    N = l - max(abs(m), abs(s)) + 1
    return N + 10
end

function _compute_normalization_constant(coefficients_params, coefficients)
    #=
    Numerically compute the integral

    \int_{0}^{pi} [nf*S(theta)]^2 sin(theta) d theta = 1

    where nf is the normalization constant needed
    =#
    
    integrand(theta) = 2*pi * _unnormalized_spin_weighted_spheroidal_harmonic(coefficients_params, coefficients, theta, 0.0)^2 * sin(theta)
    norm_sq = quadgk(theta -> integrand(theta), 0, pi)[1]
    return sqrt(norm_sq)
end

@doc raw"""
    spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c; N::Int=-1)

Construct the spectral decomposition of this spin-weighted spheroidal harmonic of 
spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$) 
using `N` spin-weighted *spherical* harmonics.

Note that the default value for `N=-1` indicates that a suitable value of `N`
will be determined automatically.

Return a SpinWeightedSpheroidalHarmonicFunction object that can be evaluated at any point.
"""
function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c; N::Int=-1)
    if N == -1
        N = _determine_matrix_size_N(s, l, m)
    end
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    coefficients = spectral_coefficients(c, s, l, m, N)
    normalization = _compute_normalization_constant(coefficients_params, coefficients)
    lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, c, N=N) # Not the most efficient way to do this, but it works

    return SpinWeightedSpheroidalHarmonicFunction(coefficients_params, coefficients, normalization, lambda)
end

# The power of multiple dispatch
@doc raw"""
    SpinWeightedSpheroidalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spheroidal harmonic at the point `(theta, phi)`. 
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSpheroidalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    _unnormalized_spin_weighted_spheroidal_harmonic(swsh_func.params, swsh_func.coeffs, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative) / swsh_func.normalization_const
end

@doc raw"""
    spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int)

Construct the spin-weighted spherical harmonic of 
spin weight `s`, harmonic index `l`, and azimuthal index `m`.

Return a SpinWeightedSphericalHarmonicFunction object that can be evaluated at any point.
"""
function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int)
    return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m))
end

@doc raw"""
    SpinWeightedSphericalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the exact value of the spin-weighted spherical harmonic at the point `(theta, phi)`.
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSphericalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    _nth_derivative_spherical_harmonic(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
end

@doc raw"""
    spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=-1)

Compute the eigenvalue of the spin-weighted spheroidal harmonic
with spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$).

The optional argument `N` specifies the number of terms to use in the spectral decomposition.
The default value is `N=-1`, which indicates that a suitable value of `N` will be determined automatically.

This function is simply a wrapper to `Teukolsky_lambda_const` for backward compatbility.
"""
function spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=-1)
    if N == -1
        N = _determine_matrix_size_N(s, l, m)
    end
    Teukolsky_lambda_const(c, s, l, m, N)
end

@doc raw"""
    spin_weighted_spherical_eigenvalue(s::Int, l::Int, m::Int=0)

Compute the eigenvalue of the spin-weighted spherical harmonic
with spin weight `s`, harmonic index `l`, and azimuthal index `m`
(but the eigenvalue is independent of `m`).

"""
function spin_weighted_spherical_eigenvalue(s::Int, l::Int, m::Int=0)
    # Eigenvalue for the Schwarzschild case does not depend on m
    Teukolsky_lambda_const(0, s, l, m)
end

end
