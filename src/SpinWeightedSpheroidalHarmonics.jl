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

function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi, derivative::Int=0, N::Int=10)
    local coefficients

    # perform spectral decomposition if inputs differ from cached params
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    if coefficients_params == _cached_coefficients_params
        coefficients = _cached_coefficients
    else
        coefficients = SpectralDecompositionCoefficients(coefficients_params)
    end

    # compute the spin-weighted spherical harmonics needed
    spherical_harmonics = []
    spherical_harmonics_l = construct_all_l_in_matrix(s, m, N)

    for spherical_harmonic_l in spherical_harmonics_l
        push!(spherical_harmonics, _nth_derivative_spherical_harmonic(s, spherical_harmonic_l, m, derivative, theta, phi))
    end

    # perform a dot product of the two vectors
    dot(coefficients, spherical_harmonics)
end

end
