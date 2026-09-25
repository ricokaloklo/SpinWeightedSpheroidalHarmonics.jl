module SpinWeightedSpheroidalHarmonics

using LinearAlgebra

include("harmonic.jl")
include("spectral.jl")
include("continuation.jl")

export spin_weighted_spheroidal_harmonic, spin_weighted_spherical_harmonic, spin_weighted_spheroidal_eigenvalue, spin_weighted_spherical_eigenvalue # Expose these functions to the user
export Teukolsky_lambda_const # For backward compatbility
export AngularEigenpair, AngularPathResult
export AngularContinuationError, AngularCache, DEFAULT_ANGULAR_CACHE
export clear_angular_cache!, continue_angular_mode, track_angular_mode
export continue_angular_lateral_pair, mirror_angular_parameters
export angular_mirror_residual, angular_precision_certificate
export angular_observables

_TOLERANCE = 1e-16 # Spherical harmonics smaller than this will be ignored in the spectral decomposition

function _format_method_name(method)
    normalized = lowercase(strip(String(method)))
    normalized in ("auto", "direct", "chebyshev", "jacobi") && return normalized
    error("Does not understand method $method. Supported values are auto, direct, chebyshev, jacobi (case-insensitive).")
end

struct SpinWeightedSphericalHarmonicFunction
    s::Int
    l::Int
    m::Int
    lambda
    method::Symbol
    chebyshev_solution
end

# Implement pretty printing for SpinWeightedSphericalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSphericalHarmonicFunction)
    print(io, "SpinWeightedSphericalHarmonicFunction(s = $(swsh_func.s), l = $(swsh_func.l), m = $(swsh_func.m), method = $(swsh_func.method), lambda = $(swsh_func.lambda))")
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
    spherical_harmonics_l::Vector{Union{SpinWeightedSphericalHarmonicFunction, Nothing}}
    normalization_const
    lambda
    method::Symbol
    chebyshev_solution
end

# Implement pretty printing for SpinWeightedSpheroidalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSpheroidalHarmonicFunction)
    print(io, "SpinWeightedSpheroidalHarmonicFunction(s = $(swsh_func.params.s), l = $(swsh_func.params.l), m = $(swsh_func.params.m), c = $(swsh_func.params.c), method = $(swsh_func.method), lambda = $(swsh_func.lambda))")
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
    return output
end

function _spheroidal_boundary_values(coefficients_params, coefficients)
    s = coefficients_params.s
    m = coefficients_params.m
    l_list = construct_all_l_in_matrix(s, m, coefficients_params.N)

    S0 = zero(eltype(coefficients))
    Spi2 = zero(eltype(coefficients))
    Spi = zero(eltype(coefficients))

    for idx in eachindex(l_list)
        if abs(coefficients[idx]) < _TOLERANCE
            continue
        end

        l = l_list[idx]
        Y0 = m == -s ? (-1)^s * sqrt((2 * l + 1) / (4π)) : 0.0
        Ypi2 = Float64(spin_weighted_spherical_harmonic_at_pi_over_2(s, l, m))
        Ypi = m == s ? (-1)^l * sqrt((2 * l + 1) / (4π)) : 0.0

        S0 += coefficients[idx] * Y0
        Spi2 += coefficients[idx] * Ypi2
        Spi += coefficients[idx] * Ypi
    end

    return S0, Spi2, Spi
end

@doc raw"""
    spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c; N::Int=-1, method="auto", backend="auto")

Construct the spectral decomposition of this spin-weighted spheroidal harmonic of 
spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$) 
using `N` spin-weighted *spherical* harmonics.

Note that the default value for `N=-1` indicates that a suitable value of `N`
will be determined automatically.

Return a SpinWeightedSpheroidalHarmonicFunction object that can be evaluated at any point.

The `method` options are `"auto"` (default), `"direct"`, `"jacobi"`, and
`"chebyshev"`. Method names are case-insensitive.

The `backend` argument controls how the spectral decomposition is solved:
- `"auto"` (default): `"banded"` for real `c`; for complex `c`, follow the mode along the path with
  Newton steps, falling back to a full dense eigendecomposition when a step is rejected,
- `"banded"`: real `c` only. Solve for the requested eigenpair with a banded solver, increasing `N`
  until it converges, which avoids a full dense eigendecomposition,
- `"dense"`: a full dense eigendecomposition. For real `c` the mode is picked by its position among
  the sorted eigenvalues; for complex `c` it is followed along the path, keeping the eigenvector with
  the largest overlap with the previous step.
Backend names are case-insensitive.

For complex `c`, "mode `l`" is the mode that connects continuously to the spherical harmonic `l`
along the straight path from `c = 0`. To follow a different path, use `track_angular_mode` and pass
a returned eigenpair to `spin_weighted_spheroidal_harmonic`.
"""
function spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c;
        N::Int=-1, method="auto", backend="auto")
    selected_backend = _resolve_spectral_backend(backend, c)
    if c isa Complex && !iszero(imag(c)) && selected_backend != :dense
        pair = continue_angular_mode(
            s, l, m, c; truncation_order=_complex_truncation_order(s, l, m, N))
        return spin_weighted_spheroidal_harmonic(pair; method)
    end
    adaptive_pair = nothing
    if N == -1 && selected_backend != :dense && isreal(c)
        adaptive_pair = _adaptive_real_eigenpair(real(c), s, l, m)
    end
    method = _format_method_name(method)
    angular_sep, coefficients = if adaptive_pair !== nothing
        shift = muladd(Float64(real(c)), Float64(real(c)),
            -2m * Float64(real(c)))
        adaptive_pair.lambda - shift, adaptive_pair.coefficients
    else
        # With N = -1, the backend picks the size, so N is read off the result
        _spectral_decomposition(c, s, l, m, N; backend=selected_backend)
    end
    N = length(coefficients)
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    normalization = 1 # already satisfied the normalization cond. \int_{0}^{pi} [nf*S(theta)]^2 sin(theta) d theta = 1
    lambda = angular_sep + c^2 - 2*m*c

    l_list = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)
    if method == "chebyshev"
        S0, Spi2, Spi = _spheroidal_boundary_values(coefficients_params, coefficients)
        chebyshev_solution = Fun(_solve_spheroidal_harmonic_chebyshev(s, m, c, lambda, S0, Spi2, Spi), 0..π)
        spherical_harmonics_l = Vector{Union{SpinWeightedSphericalHarmonicFunction, Nothing}}(undef, length(l_list))
        fill!(spherical_harmonics_l, nothing)
        return SpinWeightedSpheroidalHarmonicFunction(coefficients_params, coefficients, spherical_harmonics_l, normalization, lambda, :chebyshev, chebyshev_solution)
    else
        spherical_harmonics_l = [ abs(coefficients[n]) >= _TOLERANCE ? spin_weighted_spherical_harmonic(s, l_list[n], m; method=method) : nothing for n in eachindex(l_list) ]
        return SpinWeightedSpheroidalHarmonicFunction(coefficients_params, coefficients, spherical_harmonics_l, normalization, lambda, :spectral, nothing)
    end
end

@doc raw"""
    SpinWeightedSpheroidalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spheroidal harmonic at the point `(theta, phi)`. 
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSpheroidalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    if swsh_func.method == :spectral
        return _unnormalized_spin_weighted_spheroidal_harmonic(swsh_func.params, swsh_func.coeffs, swsh_func.spherical_harmonics_l, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative) / swsh_func.normalization_const
    elseif swsh_func.method == :chebyshev
        return _nth_derivative_spheroidal_harmonic_chebyshev(swsh_func.chebyshev_solution, swsh_func.params.m, theta_derivative, phi_derivative, theta, phi) / swsh_func.normalization_const
    end
end

@doc raw"""
    spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto")

Construct the spin-weighted spherical harmonic of 
spin weight `s`, harmonic index `l`, and azimuthal index `m`.

Return a SpinWeightedSphericalHarmonicFunction object that can be evaluated at any point.

By default, the method to compute the value is chosen automatically based on the value of the harmonic index `l`.
When `l < 30`, the direct evaluation method (`method="direct"`) is used where we evaluate the exact analytical solution as shown in Eq. (A8).
When `l >= 30`, we use a Jacobi polynomial recurrence (`method="jacobi"`) for better high-`l` stability and speed (see https://arxiv.org/abs/2208.03691).
The Chebyshev pseudo-spectral method (`method="chebyshev"`) remains available as an explicit option.
Method names are case-insensitive.
"""
function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto")
    method = _format_method_name(method)
    if method == "auto"
        method = l >= 30 ? "jacobi" : "direct"
    end

    if method == "direct"
        return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m), :direct, nothing)
    elseif method == "chebyshev"
        chebyshev_soln = _solve_spherical_harmonic_chebyshev(s, l, m)
        return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m), :chebyshev, chebyshev_soln)
    elseif method == "jacobi"
        return SpinWeightedSphericalHarmonicFunction(s, l, m, spin_weighted_spherical_eigenvalue(s, l, m), :jacobi, nothing)
    end
end

@doc raw"""
    SpinWeightedSphericalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spherical harmonic at the point `(theta, phi)`.
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSphericalHarmonicFunction)(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) = begin
    if swsh_func.method == :direct
        return _nth_derivative_spherical_harmonic_direct_eval(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    elseif swsh_func.method == :chebyshev
        return _nth_derivative_spheroidal_harmonic_chebyshev(swsh_func.chebyshev_solution, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    elseif swsh_func.method == :jacobi
        return _nth_derivative_spherical_harmonic_jacobi(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    else
        error("Unknown method $(swsh_func.method) for evaluating spin-weighted spherical harmonic.")
    end
end

@doc raw"""
    spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=-1)

Compute the eigenvalue of the spin-weighted spheroidal harmonic
with spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$).

The optional argument `N` specifies the number of terms to use in the spectral decomposition.
The default value is `N=-1`, which indicates that a suitable value of `N` will be determined automatically.

This function is simply a wrapper to `Teukolsky_lambda_const` for backward compatibility.
"""
function spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=-1)
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
