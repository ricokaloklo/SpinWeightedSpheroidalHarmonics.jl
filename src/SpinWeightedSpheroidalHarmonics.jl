module SpinWeightedSpheroidalHarmonics

using LinearAlgebra

export spin_weighted_spheroidal_harmonic

include("harmonic.jl")
include("spectral.jl")

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

_cached_coefficients_params = SpectralDecompositionInputParams(-2, 2, 2, 0.5+0.1im, 10)
_cached_coefficients = SpectralDecompositionCoefficients(_cached_coefficients_params)

function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi, theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)
    local coefficients

    # Perform spectral decomposition if inputs differ from cached params
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    if coefficients_params == _cached_coefficients_params
        coefficients = _cached_coefficients
    else
        coefficients = SpectralDecompositionCoefficients(coefficients_params)
    end

    # Compute the spin-weighted spherical harmonics needed
    output = 0.0
    spherical_harmonics_l = construct_all_l_in_matrix(s, m, N)

    for (idx, spherical_harmonic_l) in enumerate(spherical_harmonics_l)
        output += coefficients[idx] * _nth_derivative_spherical_harmonic(s, spherical_harmonic_l, m, theta_derivative, phi_derivative, theta, phi)
    end
    output
end

end
