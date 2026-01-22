module SpinWeightedSpheroidalHarmonics

using LinearAlgebra

include("harmonic.jl")
include("spectral.jl")

export spin_weighted_spheroidal_harmonic, spin_weighted_spherical_harmonic, spin_weighted_spheroidal_eigenvalue, spin_weighted_spherical_eigenvalue # Expose these functions to the user
export Teukolsky_lambda_const # For backward compatbility

_TOLERANCE = 1e-16 # Spherical harmonics smaller than this will be ignored in the spectral decomposition

struct SpinWeightedSphericalHarmonicFunction
    s::Int
    l::Int
    m::Int
    lambda
    chebyshev_solution
end

# Implement pretty printing for SpinWeightedSphericalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSphericalHarmonicFunction)
    print(io, "SpinWeightedSphericalHarmonicFunction(s = $(swsh_func.s), l = $(swsh_func.l), m = $(swsh_func.m), lambda = $(swsh_func.lambda))")
end

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
    spherical_harmonics_l::Vector{SpinWeightedSphericalHarmonicFunction}
    normalization_const
    lambda
end

# Implement pretty printing for SpinWeightedSpheroidalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSpheroidalHarmonicFunction)
    print(io, "SpinWeightedSpheroidalHarmonicFunction(s = $(swsh_func.params.s), l = $(swsh_func.params.l), m = $(swsh_func.params.m), c = $(swsh_func.params.c), lambda = $(swsh_func.lambda))")
end

function _unnormalized_spin_weighted_spheroidal_harmonic(coefficients_params, coefficients, spherical_harmonics_l, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)
    # Compute the spin-weighted spherical harmonics needed
    output = 0.0
    l_list = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)

    for (idx, _) in enumerate(l_list)
        if abs(coefficients[idx]) < _TOLERANCE
            continue
        end
        output += coefficients[idx] * spherical_harmonics_l[idx](theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative)
    end
    output
end

@doc raw"""
    spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c; N::Int=-1, method="auto")

Construct the spectral decomposition of this spin-weighted spheroidal harmonic of 
spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$) 
using `N` spin-weighted *spherical* harmonics.

Note that the default value for `N=-1` indicates that a suitable value of `N`
will be determined automatically.

Return a SpinWeightedSpheroidalHarmonicFunction object that can be evaluated at any point.

By default, the method to compute the value for the underlying spin-weighted spherical harmonics 
is chosen automatically based on the value of the harmonic index `l`.
When `l < 30`, the direct evaluation method (`method="direct"`) is used where we evaluate the exact analytical solution as shown in Eq. (A8).
However, the prefactor in each term of the sum can be very large and thus cause overflow, while the sum itself is finite (of order 1 actually).
Therefore, when `l >= 30`, the Chebyshev pseudo-spectral method (`method="chebyshev"`) is used instead.
"""
function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c; N::Int=-1, method="auto")
    if N == -1
        N = _determine_matrix_size_N(s, l, m)
    end
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    coefficients = spectral_coefficients(c, s, l, m, N)
    l_list = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)
    spherical_harmonics_l = [
        spin_weighted_spherical_harmonic(s, l, m; method=method) for l in l_list
    ]
    normalization = 1 # already satisfied the normalization cond. \int_{0}^{pi} [nf*S(theta)]^2 sin(theta) d theta = 1
    lambda = spin_weighted_spheroidal_eigenvalue(s, l, m, c, N=N) # Not the most efficient way to do this, but it works

    return SpinWeightedSpheroidalHarmonicFunction(coefficients_params, coefficients, spherical_harmonics_l, normalization, lambda)
end

# The power of multiple dispatch
@doc raw"""
    SpinWeightedSpheroidalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spheroidal harmonic at the point `(theta, phi)`. 
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSpheroidalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    _unnormalized_spin_weighted_spheroidal_harmonic(swsh_func.params, swsh_func.coeffs, swsh_func.spherical_harmonics_l, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative) / swsh_func.normalization_const
end

@doc raw"""
    spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto")

Construct the spin-weighted spherical harmonic of 
spin weight `s`, harmonic index `l`, and azimuthal index `m`.

Return a SpinWeightedSphericalHarmonicFunction object that can be evaluated at any point.

By default, the method to compute the value is chosen automatically based on the value of the harmonic index `l`.
When `l < 30`, the direct evaluation method (`method="direct"`) is used where we evaluate the exact analytical solution as shown in Eq. (A8).
However, the prefactor in each term of the sum can be very large and thus cause overflow, while the sum itself is finite (of order 1 actually).
Therefore, when `l >= 30`, the Chebyshev pseudo-spectral method (`method="chebyshev"`) is used instead.
"""
function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto")
    if method == "auto"
        _method = l >= 30 ? "chebyshev" : "direct"
    elseif method == "direct" || method == "chebyshev"
        _method = method
    else
        error("Does not understand method $method")
    end

    if _method == "direct"
        return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m), nothing)
    else
        chebyshev_soln = _solve_spherical_harmonic_chebyshev(s, l, m)
        return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m), chebyshev_soln)
    end
end

@doc raw"""
    SpinWeightedSphericalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spherical harmonic at the point `(theta, phi)`.
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSphericalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    if isnothing(swsh_func.chebyshev_solution)
        return _nth_derivative_spherical_harmonic_direct_eval(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    else
        return _nth_derivative_spherical_harmonic_chebyshev(swsh_func.chebyshev_solution, swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    end
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
