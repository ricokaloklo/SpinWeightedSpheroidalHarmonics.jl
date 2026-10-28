module SpinWeightedSpheroidalHarmonics

using LinearAlgebra

include("Precision.jl")
include("harmonic.jl")
include("spectral.jl")
include("continuation.jl")
include("leaver.jl")

export spin_weighted_spheroidal_harmonic, spin_weighted_spherical_harmonic, spin_weighted_spheroidal_eigenvalue, spin_weighted_spherical_eigenvalue
export Teukolsky_lambda_const # For backward compatibility
export AngularEigenpair, AngularPathResult
export AngularContinuationError, AngularCache, DEFAULT_ANGULAR_CACHE
export clear_angular_cache!, continue_angular_mode, track_angular_mode
export continue_angular_lateral_pair, mirror_angular_parameters
export angular_mirror_residual, angular_precision_certificate
export angular_observables
export value_and_theta_derivative

# Spherical harmonics whose coefficient is below this in magnitude are left out when the spheroidal
# harmonic is built and evaluated
_spherical_coefficient_tol(::Type{T}) where {T} = _tol(T, 1e-16)

# c (or another input) in the working type T
_to_float(::Type{T}, x::Real) where {T} = T(x)
_to_float(::Type{T}, x::Complex) where {T} = Complex{T}(x)
_to_float(::Type, ::Nothing) = nothing

function _format_method_name(method)
    normalized = lowercase(strip(String(method)))
    normalized in ("auto", "direct", "chebyshev", "jacobi", "leaver") && return normalized
    error("Does not understand method $method. Supported values are auto, direct, chebyshev, jacobi, leaver (case-insensitive).")
end

function _format_normalization(normalization)
    normalized = Symbol(lowercase(strip(String(normalization))))
    normalized in (:modulus, :analytic) && return normalized
    throw(ArgumentError("Does not understand normalization $normalization. Supported values are :modulus and :analytic."))
end

#=
The coefficients of a harmonic rescaled to the analytic normalization: since sYlm(theta, 0) is real,
\int S^2 dΩ (with S(theta, 0), not conj(S) S) is the sum of the squared coefficients, which is made 1.
Of the two square roots, the one that keeps the overlap with sYlm in the right half-plane is taken,
continuous with real c, where the coefficients are real and are returned as they are.
=#
function _analytically_normalized(coefficients::AbstractVector, s::Int, l::Int, m::Int)
    all(isreal, coefficients) && return coefficients
    scale = sqrt(sum(x -> x^2, coefficients))
    real(coefficients[l - max(abs(s), abs(m)) + 1] / scale) < 0 && (scale = -scale)
    return coefficients ./ scale
end

# T is the floating-point type the harmonic was built in; evaluation runs in the wider of T and the
# type of the angles, at the BigFloat precision the harmonic was built with (0 when T carries it)
struct SpinWeightedSphericalHarmonicFunction{T<:AbstractFloat}
    s::Int
    l::Int
    m::Int
    lambda
    method::Symbol
    chebyshev_solution
    jacobi_theta_terms::Union{Nothing, NTuple{2, Vector{Pair{NTuple{5, Int}, Float64}}}} # first and second theta-derivative term lists (Jacobi method), built once with the harmonic; the coefficients are exact
    direct_terms::Union{Nothing, DirectSumTerms{T}} # the prepared direct sum (direct method)
    precision::Int
end

# Implement pretty printing for SpinWeightedSphericalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSphericalHarmonicFunction)
    print(io, "SpinWeightedSphericalHarmonicFunction(s = $(swsh_func.s), l = $(swsh_func.l), m = $(swsh_func.m), method = $(swsh_func.method), lambda = $(swsh_func.lambda))")
end

struct SpectralDecompositionInputParams
    s::Int
    l # Not restricted to an Int, so that a future solver for non-integer degree can use it
    m::Int
    c
    N::Int
end

struct SpinWeightedSpheroidalHarmonicFunction{T<:AbstractFloat}
    params::SpectralDecompositionInputParams
    coeffs
    spherical_harmonics_l::Vector{Union{SpinWeightedSphericalHarmonicFunction, Nothing}}
    normalization_const
    lambda
    method::Symbol
    chebyshev_solution
    leaver_solution
    precision::Int
end

#=
Evaluate f at the angles converted to the evaluation type: the wider of the type T the harmonic was
built in and the type of the angles, at the precision the harmonic was built with (or that of the
angles, if wider), whatever precision the caller has set.
=#
function _evaluate_at(f, ::Type{T}, built_precision::Int, theta, phi) where {T}
    R = _float_type(zero(T), theta, phi)
    return _with_precision(R, max(built_precision, _input_precision(theta, phi))) do
        f(R(theta), R(phi))
    end
end

# Implement pretty printing for SpinWeightedSpheroidalHarmonicFunction
function Base.show(io::IO, ::MIME"text/plain", swsh_func::SpinWeightedSpheroidalHarmonicFunction)
    print(io, "SpinWeightedSpheroidalHarmonicFunction(s = $(swsh_func.params.s), l = $(swsh_func.params.l), m = $(swsh_func.params.m), c = $(swsh_func.params.c), method = $(swsh_func.method), lambda = $(swsh_func.lambda))")
end

function _unnormalized_spin_weighted_spheroidal_harmonic(::Type{T}, coefficients_params, coefficients, spherical_harmonics_l, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) where {T}
    # Sum the stored spherical harmonics weighted by their coefficients
    output = zero(theta)
    l_list = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)

    for (idx, _) in enumerate(l_list)
        # the harmonics whose coefficients are below the cutoff were left out when the harmonic was built
        if spherical_harmonics_l[idx] === nothing
            continue
        end
        output += coefficients[idx] * spherical_harmonics_l[idx](theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative)
    end
    return output
end

# S and dS/dtheta at theta = π/2, from the spectral decomposition
function _spheroidal_equator_values(::Type{T}, coefficients_params, coefficients) where {T}
    s = coefficients_params.s
    m = coefficients_params.m
    l_list = construct_all_l_in_matrix(s, m, coefficients_params.N)

    S_half = zero(eltype(coefficients))
    dS_half = zero(eltype(coefficients))

    for idx in eachindex(l_list)
        if abs(coefficients[idx]) < _spherical_coefficient_tol(T)
            continue
        end

        l = l_list[idx]
        # Evaluate sYlm(\pi/2, 0) using a numerically stable method
        Ypi2 = _spherical_harmonic_at_pi_over_2(T, s, l, m)
        dYpi2 = _spherical_harmonic_at_pi_over_2(T, s, l, m; theta_derivative=1)

        S_half += coefficients[idx] * Ypi2
        dS_half += coefficients[idx] * dYpi2
    end

    return S_half, dS_half
end

#=
Every method implemented here needs an integer harmonic index l >= max(|s|, |m|).
A non-integer (complex) degree needs a dedicated solver, which is planned; until
then it is rejected rather than silently rounded to the nearest integer mode.
=#
function _resolve_spheroidal_method(method, s::Int, l, m::Int)
    method = _format_method_name(method)
    l isa Integer || throw(ArgumentError(
        "The harmonic index l must be an integer, got l = $l. " *
        "A solver for non-integer (complex) degree is planned but not available yet."))
    lmin = max(abs(s), abs(m))
    l >= lmin || throw(ArgumentError(
        "l = $l is below the lowest mode l = max(|s|, |m|) = $lmin for s = $s, m = $m."))
    return method
end

@doc raw"""
    spin_weighted_spheroidal_harmonic(s::Int, l, m::Int, c; N::Int=-1, method="auto", backend="auto", branch_n=nothing, lambda=nothing, lambda0=nothing, cf_depth::Int=-1, tol=-1, max_iter::Int=80, precision=nothing, normalization=:modulus)

Construct the spin-weighted spheroidal harmonic of
spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$).

Return a SpinWeightedSpheroidalHarmonicFunction object that can be evaluated at any point.

The `method` argument controls how the harmonic is computed:
- `"auto"` (default): spectral decomposition, with the spherical harmonics evaluated by the Jacobi recurrence,
- `"direct"` or `"jacobi"`: spectral decomposition with that spherical-harmonic backend,
- `"chebyshev"`: solve the spheroidal ODE directly with Chebyshev pseudo-spectral collocation,
- `"leaver"`: Leaver's continued-fraction method.
Method names are case-insensitive.

For the spectral decomposition, `N` is the number of spin-weighted *spherical* harmonics used,
while for `"leaver"` it is the number of terms kept in the power series. In either case the
default value `N=-1` indicates that a suitable value of `N` will be determined automatically.

For every method, the harmonic index `l` must be an integer with `l ≥ max(|s|, |m|)`. A solver for
non-integer (complex) degree is planned.

The keyword arguments `lambda`, `lambda0`, `branch_n`, `cf_depth`, `tol` and `max_iter` are specific
to `"leaver"`. `lambda` supplies the eigenvalue
instead of solving for it, and `lambda0` starts the solver from another initial guess. `branch_n` (counted from `0`
for the lowest mode `l = max(|s|, |m|)`, and by default the mode `l` itself) chooses
which inversion of the continued fraction is used; all of them have the same roots, but the
inversion belonging to the mode is the best conditioned one. `cf_depth` sets the depth at which
the continued fraction is truncated, and `tol` and `max_iter` control its Newton iteration; both
`cf_depth=-1` and `tol=-1` pick a suitable value automatically, the latter following the
floating-point precision of `c`.

Every method fixes the overall phase by making the overlap with ${}_s Y_{lm}$ real and positive, so for
real and complex `c` alike the methods return the same function. (The eigenpairs returned by
`continue_angular_mode` and `track_angular_mode` keep the phase of the continuation along the path; a
harmonic built from one has its phase fixed in the same way.)

The `normalization` argument chooses between two normalizations, which coincide for real `c`:
- `:modulus` (default): ``\int |S|^2 \, d\Omega = 1``. For complex `c` the normalization constant
  then involves ``\bar{S}`` and is not an analytic function of `c`, so neither is the harmonic.
- `:analytic`: ``\int S^2 \, d\Omega = 2\pi \int_0^\pi S(\theta, 0)^2 \sin\theta \, d\theta = 1``,
  without complex conjugation. The harmonic is then the analytic continuation in `c` of the harmonic
  for real `c` (the overlap with ${}_s Y_{lm}$ is the continuation of a positive number, no longer
  real), which is what quasinormal-mode residues, with complex ``c = a\omega``, need.
The eigenvalue does not depend on the normalization.

Two things limit `"leaver"` in a given floating-point precision, and it raises an error rather
than return a result it cannot resolve in either case. Roundoff in the continued fraction limits
how sharply the eigenvalue is located, increasingly so as `|c|` grows, and cancellation between
the (large, alternating) power-series coefficients limits the accuracy of the harmonic itself,
increasingly so as `l` grows; the θ-derivatives, which are taken exactly from the series, share
the accuracy of the values. Both are cured by a wider floating-point precision, which the whole
method supports — pass, say, `big(c)` or `precision=256`. The spectral decomposition suffers from neither.

The `backend` argument controls how the spectral decomposition is solved (it does not apply to
`"leaver"`, which does not use a spectral decomposition):
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

Every method computes in the floating-point type of `c` (at least Float64; for `"leaver"` also of
`lambda` and `lambda0`). With `precision = p`, `c` is converted to BigFloat at `p` bits and the
harmonic is computed and evaluated at that precision. A harmonic built in BigFloat remembers its
precision: it evaluates at that precision whatever precision is in effect when it is called.
"""
function spin_weighted_spheroidal_harmonic(s::Int, l, m::Int, c;
        N::Int=-1, method="auto", backend="auto",
        branch_n::Union{Nothing, Int}=nothing, lambda=nothing, lambda0=nothing,
        cf_depth::Int=-1, tol=-1, max_iter::Int=80, precision::Union{Nothing, Integer}=nothing,
        normalization=:modulus)

    method = _resolve_spheroidal_method(method, s, l, m)
    normalization = _format_normalization(normalization)
    _check_precision_keyword(precision)
    c, lambda, lambda0 = (_with_precision_keyword(x, precision) for x in (c, lambda, lambda0))
    T = method == "leaver" ? _float_type(c, filter(!isnothing, (lambda, lambda0))...) : _float_type(c)
    bits = precision === nothing ? _input_precision(c, lambda, lambda0) : precision
    return _with_precision(T, bits) do
        _spin_weighted_spheroidal_harmonic(T, bits, s, l, m, _to_float(T, c), method, N, backend,
            branch_n, lambda, lambda0, cf_depth, tol, max_iter, normalization)
    end
end

function _spin_weighted_spheroidal_harmonic(::Type{T}, bits::Int, s::Int, l, m::Int, c, method, N, backend,
        branch_n, lambda, lambda0, cf_depth, tol, max_iter, normalization_convention) where {T}
    if method == "leaver"
        _spectral_backend(backend) == :auto || throw(ArgumentError(
            "backend=$backend selects a spectral-decomposition solver, which method=\"leaver\" does not use."))
        branch = isnothing(branch_n) ? _leaver_branch_index(s, m, l) : branch_n
        lambda = isnothing(lambda) ? _leaver_eigenvalue(s, l, m, c; branch_n=branch, lambda0=lambda0, cf_depth=cf_depth, tol=tol, max_iter=max_iter) : lambda
        leaver_solution = _leaver_angular_solution(s, m, c, lambda; branch_n=branch, nmax=N)
        normalization = _leaver_normalization_const(leaver_solution, s, l, m;
            analytic=normalization_convention == :analytic)
        coefficients_params = SpectralDecompositionInputParams(s, l, m, c, length(leaver_solution.coeffs))
        return SpinWeightedSpheroidalHarmonicFunction{T}(coefficients_params, leaver_solution.coeffs, Vector{Union{SpinWeightedSphericalHarmonicFunction, Nothing}}(), normalization, lambda, :leaver, nothing, leaver_solution, bits)
    end

    selected_backend = _resolve_spectral_backend(backend, c)
    if c isa Complex && !iszero(imag(c)) && selected_backend != :dense
        pair = continue_angular_mode(
            s, l, m, c; truncation_order=_complex_truncation_order(s, l, m, N, c))
        return spin_weighted_spheroidal_harmonic(pair; method, normalization=normalization_convention)
    end
    adaptive_pair = nothing
    if N == -1 && selected_backend != :dense && isreal(c)
        adaptive_pair = _adaptive_real_eigenpair(real(c), s, l, m)
    end
    angular_sep, coefficients = if adaptive_pair !== nothing
        shift = muladd(T(real(c)), T(real(c)),
            -2m * T(real(c)))
        adaptive_pair.lambda - shift, adaptive_pair.coefficients
    else
        # With N = -1, the backend picks the size, so N is read off the result
        _spectral_decomposition(c, s, l, m, N; backend=selected_backend)
    end
    N = length(coefficients)
    if normalization_convention == :analytic
        coefficients = _analytically_normalized(coefficients, s, l, m)
    end
    coefficients_params = SpectralDecompositionInputParams(s, l, m, c, N)
    normalization = 1 # the coefficients have unit norm, so that \int_{0}^{pi} |S(theta)|^2 sin(theta) d theta = 1/(2pi) (or S^2, with the analytic normalization)
    lambda = angular_sep + c^2 - 2*m*c

    l_list = construct_all_l_in_matrix(coefficients_params.s, coefficients_params.m, coefficients_params.N)
    if method == "chebyshev"
        S_half, dS_half = _spheroidal_equator_values(T, coefficients_params, coefficients)
        chebyshev_solution = _solve_spheroidal_harmonic_chebyshev(T, s, l, m, c, angular_sep, S_half, dS_half)
        spherical_harmonics_l = Vector{Union{SpinWeightedSphericalHarmonicFunction, Nothing}}(undef, length(l_list))
        fill!(spherical_harmonics_l, nothing)
        return SpinWeightedSpheroidalHarmonicFunction{T}(coefficients_params, coefficients, spherical_harmonics_l, normalization, lambda, :chebyshev, chebyshev_solution, nothing, bits)
    else
        spherical_harmonics_l = [ abs(coefficients[n]) >= _spherical_coefficient_tol(T) ? spin_weighted_spherical_harmonic(s, l_list[n], m; method, precision=_precision_keyword(bits)) : nothing for n in eachindex(l_list) ]
        return SpinWeightedSpheroidalHarmonicFunction{T}(coefficients_params, coefficients, spherical_harmonics_l, normalization, lambda, :spectral, nothing, nothing, bits)
    end
end

@doc raw"""
    SpinWeightedSpheroidalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spheroidal harmonic at the point `(theta, phi)`. 
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSpheroidalHarmonicFunction{T})(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) where {T} =
    _evaluate_at(T, swsh_func.precision, theta, phi) do theta, phi
    if swsh_func.method == :spectral
        return _unnormalized_spin_weighted_spheroidal_harmonic(T, swsh_func.params, swsh_func.coeffs, swsh_func.spherical_harmonics_l, theta, phi; theta_derivative=theta_derivative, phi_derivative=phi_derivative) / swsh_func.normalization_const
    elseif swsh_func.method == :chebyshev
        return _nth_derivative_spheroidal_harmonic_chebyshev(swsh_func.chebyshev_solution, swsh_func.params.m, theta_derivative, phi_derivative, theta, phi) / swsh_func.normalization_const
    elseif swsh_func.method == :leaver
        return _nth_derivative_spheroidal_harmonic_leaver(swsh_func.leaver_solution, swsh_func.params.m, theta_derivative, phi_derivative, theta, phi) / swsh_func.normalization_const
    else
        error("Unknown method $(swsh_func.method) for evaluating spin-weighted spheroidal harmonic.")
    end
end

@doc raw"""
    spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto", precision=nothing)

Construct the spin-weighted spherical harmonic of 
spin weight `s`, harmonic index `l`, and azimuthal index `m`.

Return a SpinWeightedSphericalHarmonicFunction object that can be evaluated at any point.

By default (`method="auto"`) the harmonic is evaluated with a Jacobi polynomial recurrence (`method="jacobi"`, see
https://arxiv.org/abs/2208.03691), which is stable and fast for every `l`, near the poles too. The direct evaluation of the
exact analytical sum, Eq. (A8) (`method="direct"`), can be selected explicitly; it loses digits to cancellation as `l` grows.
The Chebyshev pseudo-spectral method (`method="chebyshev"`) can be selected explicitly.
Method names are case-insensitive.

The harmonic is evaluated in the floating-point type of the angles (at least Float64). With
`precision = p` it is built and evaluated in BigFloat at `p` bits, whatever precision is in effect
when it is called.
"""
function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int; method="auto",
        precision::Union{Nothing, Integer}=nothing)
    method = _format_method_name(method)
    if method == "auto"
        method = "jacobi"
    end
    _check_precision_keyword(precision)
    T = precision === nothing ? Float64 : BigFloat
    bits = precision === nothing ? 0 : Int(precision)
    lambda = spin_weighted_spherical_eigenvalue(s, l, m)

    if method == "direct"
        direct_terms = _with_precision(() -> DirectSumTerms(T, s, l, m), T, bits)
        return SpinWeightedSphericalHarmonicFunction{T}(s, l, m, lambda, :direct, nothing, nothing, direct_terms, bits)
    elseif method == "chebyshev"
        chebyshev_soln = _with_precision(() -> _solve_spherical_harmonic_chebyshev(T, s, l, m), T, bits)
        return SpinWeightedSphericalHarmonicFunction{T}(s, l, m, lambda, :chebyshev, chebyshev_soln, nothing, nothing, bits)
    elseif method == "jacobi"
        return SpinWeightedSphericalHarmonicFunction{T}(s, l, m, lambda, :jacobi, nothing,
            _jacobi_theta_derivative_terms(s, l, m), nothing, bits)
    else
        error("Method $method is not available for spin-weighted spherical harmonics. Supported values are auto, direct, chebyshev, jacobi (case-insensitive).")
    end
end

@doc raw"""
    SpinWeightedSphericalHarmonicFunction(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the value of the spin-weighted spherical harmonic at the point `(theta, phi)`.
Additionally compute the `theta_derivative`-th derivative with respect to `theta` and the `phi_derivative`-th derivative with respect to `phi` exactly.
"""
(swsh_func::SpinWeightedSphericalHarmonicFunction{T})(theta, phi; theta_derivative::Int=0, phi_derivative::Int=0) where {T} =
    _evaluate_at(T, swsh_func.precision, theta, phi) do theta, phi
    if swsh_func.method == :direct
        return _nth_derivative_spherical_harmonic_direct_eval(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi;
            prepared=swsh_func.direct_terms)
    elseif swsh_func.method == :chebyshev
        return _nth_derivative_spheroidal_harmonic_chebyshev(swsh_func.chebyshev_solution, swsh_func.m, theta_derivative, phi_derivative, theta, phi)
    elseif swsh_func.method == :jacobi
        return _nth_derivative_spherical_harmonic_jacobi(swsh_func.s, swsh_func.l, swsh_func.m, theta_derivative, phi_derivative, theta, phi;
            prepared_terms=swsh_func.jacobi_theta_terms)
    else
        error("Unknown method $(swsh_func.method) for evaluating spin-weighted spherical harmonic.")
    end
end

@doc raw"""
    value_and_theta_derivative(harmonic, theta, phi=0)
    value_and_theta_derivative(harmonic, thetas::AbstractVector, phi=0)

Evaluate a spin-weighted spheroidal or spherical harmonic and its first derivative with respect to
`theta` together, and return them as a tuple `(S, dS/dtheta)`; for a vector of angles, return the
tuple of the two vectors. The results are those of `harmonic(theta, phi)` and
`harmonic(theta, phi; theta_derivative=1)`, in the same floating-point type and precision.

For a harmonic built with the Jacobi recurrence (the default, `method="auto"` or `"jacobi"`, of
`spin_weighted_spheroidal_harmonic` and `spin_weighted_spherical_harmonic`) the two are computed in
one pass, which is several times faster than the two separate calls: the Jacobi polynomials of all
the spherical harmonics in the expansion come from one recurrence for the value and one for the
derivative. The evaluation is accurate up to the poles. For the other methods, and at the poles
themselves, this is the same as the two separate calls.

# Example
```julia
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, 0.35)
value, derivative = value_and_theta_derivative(S, 1.1, 0.3)
values, derivatives = value_and_theta_derivative(S, range(0, π; length=101))
```
"""
function value_and_theta_derivative(harmonic::SpinWeightedSpheroidalHarmonicFunction{T}, theta, phi=0) where {T}
    # As _evaluate_at, written out so that a Float64 evaluation does not allocate
    R = _float_type(zero(T), theta, phi)
    return _with_precision(R, max(harmonic.precision, _input_precision(theta, phi))) do
        _value_and_theta_derivative(harmonic, R(theta), R(phi))::NTuple{2, Complex{R}}
    end
end

function value_and_theta_derivative(harmonic::SpinWeightedSphericalHarmonicFunction{T}, theta, phi=0) where {T}
    # As _evaluate_at, written out so that a Float64 evaluation does not allocate
    R = _float_type(zero(T), theta, phi)
    return _with_precision(R, max(harmonic.precision, _input_precision(theta, phi))) do
        _value_and_theta_derivative(harmonic, R(theta), R(phi))::NTuple{2, Complex{R}}
    end
end

function value_and_theta_derivative(harmonic::SpinWeightedSpheroidalHarmonicFunction{T}, thetas::AbstractVector, phi=0) where {T}
    return _values_and_theta_derivatives(harmonic, T, thetas, phi)
end

function value_and_theta_derivative(harmonic::SpinWeightedSphericalHarmonicFunction{T}, thetas::AbstractVector, phi=0) where {T}
    return _values_and_theta_derivatives(harmonic, T, thetas, phi)
end

# As _evaluate_at, for many angles
function _values_and_theta_derivatives(harmonic, ::Type{T}, thetas::AbstractVector, phi) where {T}
    R = _float_type(zero(T), zero(eltype(thetas)), phi)
    return _with_precision(R, max(harmonic.precision, _input_precision(thetas, phi))) do
        values = Vector{Complex{R}}(undef, length(thetas))
        derivatives = similar(values)
        _values_and_theta_derivatives!(values, derivatives, harmonic, thetas, R(phi))
        values, derivatives
    end
end

#=
The one-pass sum is used for a harmonic evaluated with the Jacobi recurrence; the arguments it is
given have concrete types (the harmonic stores its coefficients untyped), so that the angles are
evaluated without allocating. Any other harmonic is evaluated by the two separate calls.
=#
function _one_pass_harmonics(harmonic::SpinWeightedSpheroidalHarmonicFunction)
    harmonic.method == :spectral || return nothing
    harmonics = harmonic.spherical_harmonics_l
    first_index = findfirst(!isnothing, harmonics)
    first_index === nothing && return nothing
    return harmonics[first_index].method == :jacobi ? harmonics : nothing
end

function _value_and_theta_derivative(harmonic::SpinWeightedSpheroidalHarmonicFunction{T}, theta, phi) where {T}
    harmonics = _one_pass_harmonics(harmonic)
    harmonics === nothing && return _two_calls(harmonic, theta, phi)
    params = harmonic.params
    lmin = max(abs(params.s), abs(params.m))
    coefficients, normalization = harmonic.coeffs, harmonic.normalization_const
    # The types a spectral harmonic stores, spelled out so that the call below is resolved statically
    if coefficients isa Vector{Complex{T}} && normalization isa Int
        return _one_pass_at(harmonic, params.s, params.m, lmin, lmin + params.N - 1,
            coefficients, harmonics, normalization, theta, phi)
    elseif coefficients isa Vector{Complex{T}} && normalization isa T
        return _one_pass_at(harmonic, params.s, params.m, lmin, lmin + params.N - 1,
            coefficients, harmonics, normalization, theta, phi)
    end
    return _one_pass_at(harmonic, params.s, params.m, lmin, lmin + params.N - 1,
        coefficients, harmonics, normalization, theta, phi)::NTuple{2, Complex{typeof(theta)}}
end

_value_and_theta_derivative(harmonic::SpinWeightedSphericalHarmonicFunction, theta, phi) =
    harmonic.method == :jacobi ?
        _one_pass_at(harmonic, harmonic.s, harmonic.m, harmonic.l, harmonic.l, nothing, nothing, 1, theta, phi) :
        _two_calls(harmonic, theta, phi)

function _values_and_theta_derivatives!(values, derivatives, harmonic::SpinWeightedSpheroidalHarmonicFunction, thetas, phi)
    harmonics = _one_pass_harmonics(harmonic)
    harmonics === nothing && return _two_calls!(values, derivatives, harmonic, thetas, phi)
    params = harmonic.params
    lmin = max(abs(params.s), abs(params.m))
    return _one_pass_at!(values, derivatives, harmonic, params.s, params.m, lmin, lmin + params.N - 1,
        harmonic.coeffs, harmonics, harmonic.normalization_const, thetas, phi)
end

_values_and_theta_derivatives!(values, derivatives, harmonic::SpinWeightedSphericalHarmonicFunction, thetas, phi) =
    harmonic.method == :jacobi ?
        _one_pass_at!(values, derivatives, harmonic, harmonic.s, harmonic.m, harmonic.l, harmonic.l, nothing, nothing, 1, thetas, phi) :
        _two_calls!(values, derivatives, harmonic, thetas, phi)

function _two_calls(harmonic, theta, phi)
    C = Complex{typeof(theta)}
    return convert(C, harmonic(theta, phi))::C, convert(C, harmonic(theta, phi; theta_derivative=1))::C
end

function _two_calls!(values, derivatives, harmonic, thetas, phi)
    for (i, theta) in enumerate(thetas)
        values[i], derivatives[i] = _two_calls(harmonic, real(eltype(values))(theta), phi)
    end
end

function _one_pass_at!(values, derivatives, harmonic, s, m, lmin, lmax, coefficients, harmonics, normalization, thetas, phi)
    for (i, theta) in enumerate(thetas)
        values[i], derivatives[i] = _one_pass_at(harmonic, s, m, lmin, lmax, coefficients, harmonics, normalization,
            real(eltype(values))(theta), phi)
    end
end

function _one_pass_at(harmonic, s, m, lmin, lmax, coefficients, harmonics, normalization, theta, phi)
    _theta, _phi = _fold_theta_phi(theta, phi)
    # At the poles the two evaluations take their exact pole values
    _near_pole(_theta) && return _two_calls(harmonic, theta, phi)
    weight = coefficients === nothing ? (l -> 1) :
        (l -> (i = l - lmin + 1; harmonics[i] === nothing ? nothing : coefficients[i]))
    value, derivative = _jacobi_value_and_theta_derivative(weight, s, m, lmin, lmax, _theta, _phi)
    C = Complex{typeof(theta)}
    return convert(C, value / normalization), convert(C, derivative / normalization)
end

@doc raw"""
    spin_weighted_spheroidal_eigenvalue(s::Int, l, m::Int, c; N::Int=-1, method="auto", branch_n=nothing, lambda0=nothing, cf_depth::Int=-1, tol=-1, max_iter::Int=80, precision=nothing)

Compute the eigenvalue of the spin-weighted spheroidal harmonic
with spin weight `s`, harmonic index `l`, azimuthal index `m`, and spheroidicity `c` ($c = a\omega$).

The optional argument `N` specifies the number of terms to use in the spectral decomposition.
The default value is `N=-1`, which indicates that a suitable value of `N` will be determined automatically.

With `method="leaver"`, the eigenvalue is instead obtained from Leaver's continued fraction. See
[`spin_weighted_spheroidal_harmonic`](@ref) for the meaning of `branch_n`, `lambda0`, `cf_depth`,
`tol` and `max_iter`, which are all specific to that method. `method="auto"` selects the spectral
decomposition. For every method, `l` must be an integer with `l ≥ max(|s|, |m|)`.

For the spectral decomposition this function is simply a wrapper to `Teukolsky_lambda_const`
for backward compatibility.

The eigenvalue is computed in the floating-point type of `c` (at least Float64), or in BigFloat at
`precision` bits.
"""
function spin_weighted_spheroidal_eigenvalue(s::Int, l, m::Int, c; N::Int=-1, method="auto",
    branch_n::Union{Nothing, Int}=nothing, lambda0=nothing, cf_depth::Int=-1, tol=-1, max_iter::Int=80,
    precision::Union{Nothing, Integer}=nothing)

    method = _resolve_spheroidal_method(method, s, l, m)
    _check_precision_keyword(precision)
    c, lambda0 = (_with_precision_keyword(x, precision) for x in (c, lambda0))
    T = method == "leaver" ? _float_type(c, filter(!isnothing, (lambda0,))...) : _float_type(c)
    bits = precision === nothing ? _input_precision(c, lambda0) : precision
    return _with_precision(T, bits) do
        c = _to_float(T, c)
        if method == "leaver"
            return _leaver_eigenvalue(s, l, m, c; branch_n=branch_n, lambda0=lambda0, cf_depth=cf_depth, tol=tol, max_iter=max_iter)
        end
        # N = -1 is passed on: the spectral solvers then choose the matrix size themselves
        # (_adaptive_real_lambda for real c, _complex_truncation_order for complex c)
        Teukolsky_lambda_const(c, s, l, m, N)
    end
end

@doc raw"""
    spin_weighted_spherical_eigenvalue(s::Int, l::Int, m::Int=0)

Compute the eigenvalue of the spin-weighted spherical harmonic
with spin weight `s`, harmonic index `l`, and azimuthal index `m`
(but the eigenvalue is independent of `m`).

"""
function spin_weighted_spherical_eigenvalue(s::Int, l::Int, m::Int=0)
    # Eigenvalue for the Schwarzschild case does not depend on m; it is an exact integer
    eigenvalue_Schwarzschild(s, l)
end


end
