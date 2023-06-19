module SpinWeightedSpheroidalHarmonics

using LinearAlgebra
using QuadGK

include("harmonic.jl")
include("spectral.jl")

export spin_weighted_spheroidal_harmonic, spin_weighted_spherical_harmonic, spin_weighted_spheroidal_eigenvalue
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

function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)
    global _cached_coefficients_params
    global _cached_coefficients
    global _cached_normalization

    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    # Perform spectral decomposition and normalization if inputs differ from cached params
    if coefficients_params != _cached_coefficients_params
        _cached_coefficients_params = coefficients_params
        _cached_coefficients = SpectralDecompositionCoefficients(_cached_coefficients_params)
        _cached_normalization = _compute_normalization_constant(_cached_coefficients_params, _cached_coefficients)
    end

    unnormalized_harmonic = _unnormalized_spin_weighted_spheroidal_harmonic(_cached_coefficients_params, _cached_coefficients, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative)
    unnormalized_harmonic/_cached_normalization
end

function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)
    spin_weighted_spheroidal_harmonic(s, l, m, 0, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative, N=N)
end

function spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=10)
    Teukolsky_lambda_const(c, s, l, m, N)
end

end
