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

struct SpectralDecompositionCoefficients
    input::SpectralDecompositionInputParams
    SpectralDecompositionCoefficients(input) = spectral_coefficients(input.c, input.s, input.l, input.m, input.N)
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

_cached_coefficients_params = SpectralDecompositionInputParams(-2, 2, 2, 0.5+0.1im, 10)
_cached_coefficients = SpectralDecompositionCoefficients(_cached_coefficients_params)
_cached_normalization = _compute_normalization_constant(_cached_coefficients_params, _cached_coefficients)

@doc raw"""
    spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)

Compute the spin-weighted spheroidal harmonic with spin weight `s`, harmonic index `l`, azimuthal index `m`,
spheroidicity `c` ($c = a\omega$), and coordinates `theta` and `phi`.

The optional arguments `theta_derivative` and `phi_derivative` specify the order of partial derivatives to take with respect to `theta` and `phi`, respectively.

The optional argument `N` specifies the number of terms to use in the spectral decomposition. The default value is `N=10`.
"""
function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)
    global _cached_coefficients_params
    global _cached_coefficients
    global _cached_normalization

    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    # Perform spectral decomposition and normalization if inputs differ from cached params
    if coefficients_params != _cached_coefficients_params
        _cached_coefficients_params = coefficients_params
        _cached_coefficients = SpectralDecompositionCoefficients(_cached_coefficients_params)
        # Note that if everything is consistent, this should return 1.0 (so we are not really re-normalizing the harmonic)
        _cached_normalization = _compute_normalization_constant(_cached_coefficients_params, _cached_coefficients)
    end

    unnormalized_harmonic = _unnormalized_spin_weighted_spheroidal_harmonic(_cached_coefficients_params, _cached_coefficients, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative)
    unnormalized_harmonic/_cached_normalization
end

@doc raw"""
    spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=10)

Compute the eigenvalue of the spin-weighted spheroidal harmonic
with spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$).

The optional argument `N` specifies the number of terms to use in the spectral decomposition. The default value is `N=10`.

This function is simply a wrapper to `Teukolsky_lambda_const` for backward compatbility.
"""
function spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=10)
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
