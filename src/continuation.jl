# Continuation of an angular mode along a path of spheroidicities, followed in Float64 (Newton steps
# with a dense eigendecomposition to fall back on) and refined in a wider working precision; the cache
# of continued modes; the precision certificate.
#=
Defaults of the keyword arguments of track_angular_mode (see its docstring), and the limit on the
eigenvalue condition number. A step changes c by at most SWSH_DEFAULT_MAX_PATH_STEP; a rejected step is
halved until it falls below SWSH_DEFAULT_MIN_PATH_STEP relative to max(1, |c|), at which point the
rejection is taken to be a collision rather than a step that is too large. SWSH_DEFAULT_OVERLAP_MIN is
the overlap with the previous eigenvector that identifies a candidate as the same mode, and
SWSH_DEFAULT_OVERLAP_MARGIN_MIN the lead over the runner-up required when two eigenvalues are close.
SWSH_DEFAULT_RESIDUAL_TOL bounds the relative residual of an accepted eigenpair.
SWSH_DEFAULT_EP_CONDITION_LIMIT bounds the eigenvalue condition number 1/|vᵀv| (see
_complex_symmetric_eigenvector_condition), which diverges at an exceptional point: a candidate above it
is rejected, and one above its square root is flagged :near_collision.
=#
const SWSH_DEFAULT_MAX_PATH_STEP = 0.5
const SWSH_DEFAULT_MIN_PATH_STEP = 1.0e-10
const SWSH_DEFAULT_OVERLAP_MIN = 0.65
const SWSH_DEFAULT_OVERLAP_MARGIN_MIN = 1.0e-8
const SWSH_DEFAULT_RESIDUAL_TOL = 5.0e-12
const SWSH_DEFAULT_EP_CONDITION_LIMIT = 1.0e8

struct AngularCacheKey{T<:AbstractFloat}
    s::Int
    l::Int
    m::Int
    c::Complex{T}
    sheet_id::Symbol
    precision_bits::Int
    truncation_order::Int
    backend::Symbol
    options::NamedTuple # the continuation keywords, which the continued mode also depends on
end

"""
    AngularEigenpair{T}

A spin-weighted spheroidal eigenpair at a complex spheroidicity, as returned by
`continue_angular_mode` and `track_angular_mode`. Pass it to `spin_weighted_spheroidal_harmonic`
to build the harmonic. `T` is the floating-point type it was computed in.

The main fields are `s`, `l`, `m`, `c`, the eigenvalue `lambda`, `angular_sep` (`lambda - c^2 + 2mc`),
the spherical-harmonic `coefficients` (unit norm, in the phase of the continuation along the path)
and `matrix_size`. The others describe how the state was obtained: `residual` (relative residual of
the eigen-equation), `previous_overlap` and `overlap_margin` (with the previous step, and how much
larger than the runner-up's), `spectral_gap` (relative distance to the nearest other eigenvalue),
`eigenvector_condition`, `sheet_id`, `truncation_order`, `precision_bits` (the significand bits of
`T`) and `status`, how the state was obtained: `:spherical_anchor` (`c = 0`), `:predictor_corrected`
(Newton step), `:continued` (dense eigendecomposition), `:banded`, `:precision_tracked` (a step taken
in the working precision where Float64 cannot resolve the mode) or `:near_collision` (accepted, but
close to another eigenvalue or to an exceptional point).
"""
struct AngularEigenpair{T<:AbstractFloat}
    s::Int
    l::Int
    m::Int
    c::Complex{T}
    sheet_id::Symbol
    angular_sep::Complex{T}
    lambda::Complex{T}
    coefficients::Vector{Complex{T}}
    matrix_size::Int
    truncation_order::Int
    precision_bits::Int
    normalization::Symbol
    residual::T
    previous_overlap::T
    overlap_margin::T
    spectral_gap::T
    eigenvector_condition::T
    phase_factor::Complex{T}
    status::Symbol
end

"""
    AngularPathResult{T}

The result of `track_angular_mode`: `states` (an `AngularEigenpair` for `c = 0`, every intermediate
step and every point of the path), the `requested_path`, its `sheet_id`, and for a closed path
`monodromy_overlap` and `lambda_closure_error`. `status` is `:open`, `:closed` or `:monodromy_failed`.
"""
struct AngularPathResult{T<:AbstractFloat}
    states::Vector{AngularEigenpair{T}}
    requested_path::Vector{Complex{T}}
    sheet_id::Symbol
    monodromy_overlap::Complex{T}
    lambda_closure_error::T
    status::Symbol
end

"""
    AngularContinuationError

Thrown when a step of the continuation cannot be accepted even at the smallest step size, which
usually means the path runs into an eigenvalue collision or an exceptional point. It records the
step (`c_from`, `c_to`) and the diagnostics of the last attempt.
"""
struct AngularContinuationError <: Exception
    message::String
    c_from::ComplexF64
    c_to::ComplexF64
    overlap::Float64
    overlap_margin::Float64
    residual::Float64
    spectral_gap::Float64
    eigenvector_condition::Float64
end

function Base.showerror(io::IO, exception::AngularContinuationError)
    print(io, exception.message,
        " c_from=", exception.c_from,
        " c_to=", exception.c_to,
        " overlap=", exception.overlap,
        " overlap_margin=", exception.overlap_margin,
        " residual=", exception.residual,
        " spectral_gap=", exception.spectral_gap,
        " eigenvector_condition=", exception.eigenvector_condition)
end

"""
    AngularCache(; max_entries=4096)

A thread-safe cache of continued modes, used by `continue_angular_mode`. When it holds
`max_entries` results, the oldest one is dropped. `DEFAULT_ANGULAR_CACHE` is the cache used by
default.
"""
mutable struct AngularCache
    values::Dict{Any,Any}
    insertion_order::Vector{Any}
    max_entries::Int
    lock::ReentrantLock
end

function AngularCache(; max_entries::Int=4096)
    max_entries > 0 || throw(ArgumentError("max_entries must be positive."))
    return AngularCache(Dict{Any,Any}(), Any[], max_entries, ReentrantLock())
end

const DEFAULT_ANGULAR_CACHE = AngularCache()

"""
    clear_angular_cache!(cache=DEFAULT_ANGULAR_CACHE)

Empty `cache` and return it.
"""
function clear_angular_cache!(cache::AngularCache=DEFAULT_ANGULAR_CACHE)
    lock(cache.lock) do
        empty!(cache.values)
        empty!(cache.insertion_order)
    end
    return cache
end

function _angular_cache_get(cache::AngularCache, key)
    return lock(cache.lock) do
        get(cache.values, key, nothing)
    end
end

function _angular_cache_put!(cache::AngularCache, key, value)
    lock(cache.lock) do
        if !haskey(cache.values, key)
            push!(cache.insertion_order, key)
        end
        cache.values[key] = value
        while length(cache.insertion_order) > cache.max_entries
            oldest = popfirst!(cache.insertion_order)
            delete!(cache.values, oldest)
        end
    end
    return value
end

@inline function _angular_matrix_size(
    s::Int,
    l::Int,
    m::Int,
    truncation_order::Int,
)
    truncation_order >= 0 ||
        throw(ArgumentError("truncation_order must be nonnegative."))
    index = l - max(abs(m), abs(s)) + 1
    index >= 1 || throw(ArgumentError("l is outside the angular basis."))
    return index + truncation_order
end

function _resize_coefficients(coefficients, size::Int)
    result = zeros(eltype(coefficients), size)
    count = min(length(coefficients), size)
    copyto!(result, 1, coefficients, 1, count)
    return result
end

function _phase_align!(coefficients, reference)
    resized_reference = _resize_coefficients(reference, length(coefficients))
    overlap = dot(resized_reference, coefficients)
    if !iszero(overlap)
        phase = overlap / abs(overlap)
        coefficients ./= phase
        return phase, abs(overlap)
    end
    pivot = argmax(abs.(coefficients))
    phase = iszero(coefficients[pivot]) ? one(eltype(coefficients)) :
        coefficients[pivot] / abs(coefficients[pivot])
    coefficients ./= phase
    return phase, zero(real(eltype(coefficients)))
end

function _relative_eigen_residual(matrix, value, coefficients)
    residual = norm(matrix * coefficients - value * coefficients)
    scale = (opnorm(matrix, Inf) + abs(value)) * norm(coefficients)
    return residual / max(scale, eps(real(eltype(matrix))))
end

# Gap to the nearest other eigenvalue, relative to the matrix norm: this is what
# decides whether the eigenvector can be resolved (see _eigenvector_gap_tol)
function _relative_spectral_gap(values, index::Int, matrix_norm)
    length(values) == 1 && return Inf
    value = values[index]
    gap = minimum(abs(values[j] - value) for j in eachindex(values)
        if j != index)
    return gap / max(matrix_norm, floatmin(typeof(matrix_norm)))
end

function _spherical_angular_pair(
    s::Int,
    l::Int,
    m::Int,
    matrix_size::Int,
    truncation_order::Int,
    sheet_id::Symbol,
)
    index = _ell_index_in_matrix(s, l, m, matrix_size)
    coefficients = zeros(ComplexF64, matrix_size)
    coefficients[index] = 1.0 + 0.0im
    value = ComplexF64(eigenvalue_Schwarzschild(s, l))
    return AngularEigenpair(
        s, l, m, 0.0 + 0.0im, sheet_id, value, value,
        coefficients, matrix_size, truncation_order, 53,
        :parallel_transport, 0.0, 1.0, 1.0, Inf, 1.0,
        1.0 + 0.0im,
        :spherical_anchor)
end

# For a complex-symmetric matrix the left eigenvector is the transpose of the right one, so 1/|vᵀv|
# (unit v) is the condition number of the eigenvalue: 1 for a real eigenvector, and unbounded at an
# exceptional point, where v becomes self-orthogonal (floored at 1/eps here)
@inline function _complex_symmetric_eigenvector_condition(coefficients)
    self_overlap = abs(sum(value * value for value in coefficients))
    return inv(max(self_overlap, eps(real(eltype(coefficients)))))
end

function _angular_candidate(
    previous::AngularEigenpair{Float64},
    c::ComplexF64,
)
    matrix = construct_spectral_matrix(
        c, previous.s, previous.m, previous.matrix_size)
    decomposition = eigen(matrix)
    reference = _resize_coefficients(
        previous.coefficients, previous.matrix_size)
    overlaps = Vector{Float64}(undef, length(decomposition.values))
    @inbounds for index in eachindex(decomposition.values)
        vector = view(decomposition.vectors, :, index)
        overlaps[index] = abs(dot(reference, vector)) /
            max(norm(reference) * norm(vector), floatmin(Float64))
    end
    permutation = sortperm(overlaps; rev=true)
    selected = first(permutation)
    coefficients = Vector{ComplexF64}(
        view(decomposition.vectors, :, selected))
    coefficients ./= norm(coefficients)
    phase, overlap = _phase_align!(coefficients, reference)
    second_overlap = length(permutation) > 1 ? overlaps[permutation[2]] : 0.0
    margin = overlap - second_overlap
    angular_sep = ComplexF64(decomposition.values[selected])
    lambda = angular_sep + c^2 - 2previous.m * c
    residual = _relative_eigen_residual(
        matrix, angular_sep, coefficients)
    gap = _relative_spectral_gap(
        decomposition.values, selected, opnorm(matrix, Inf))
    condition = _complex_symmetric_eigenvector_condition(coefficients)
    # The status is informational: _angular_candidate_accepted applies its own tests to the gap and the
    # condition number. Below a relative gap of sqrt(eps) roundoff has turned the eigenvector by about
    # eps/gap > sqrt(eps), i.e. it has lost half its digits; above sqrt(limit) the condition number is
    # halfway, in digits, to SWSH_DEFAULT_EP_CONDITION_LIMIT
    status = gap <= sqrt(eps(Float64)) ||
        condition >= sqrt(SWSH_DEFAULT_EP_CONDITION_LIMIT) ?
        :near_collision : :continued
    return AngularEigenpair(
        previous.s, previous.l, previous.m, c, previous.sheet_id,
        angular_sep, lambda, coefficients, previous.matrix_size,
        previous.truncation_order, 53, :parallel_transport,
        residual, overlap, margin, gap, condition, phase, status)
end

function _angular_predictor(
    states::Vector{AngularEigenpair{T}},
    previous::AngularEigenpair{T},
    target::Complex{T},
) where {T}
    length(states) < 2 &&
        return previous.angular_sep, copy(previous.coefficients)
    before = states[end - 1]
    denominator = previous.c - before.c
    iszero(denominator) &&
        return previous.angular_sep, copy(previous.coefficients)
    ratio = (target - previous.c) / denominator
    value = previous.angular_sep +
        ratio * (previous.angular_sep - before.angular_sep)
    before_coefficients = _resize_coefficients(
        before.coefficients, length(previous.coefficients))
    coefficients = previous.coefficients +
        ratio * (previous.coefficients - before_coefficients)
    iszero(norm(coefficients)) && (coefficients = copy(previous.coefficients))
    coefficients ./= norm(coefficients)
    _phase_align!(coefficients, previous.coefficients)
    return value, coefficients
end

function _bordered_jacobian(matrix, value, coefficients, reference)
    size = length(coefficients)
    jacobian = zeros(ComplexF64, size + 1, size + 1)
    jacobian[1:size, 1:size] .= matrix
    @inbounds for index in 1:size
        jacobian[index, index] -= value
        jacobian[index, size + 1] = -coefficients[index]
        jacobian[size + 1, index] = conj(reference[index])
    end
    return jacobian
end

function _corrected_angular_candidate(
    previous::AngularEigenpair{Float64},
    target::ComplexF64,
    value_seed::ComplexF64,
    coefficient_seed::Vector{ComplexF64};
    max_iterations::Int=12,
)
    matrix = construct_spectral_matrix(
        target, previous.s, previous.m, previous.matrix_size)
    coefficients = copy(coefficient_seed)
    reference = _resize_coefficients(
        previous.coefficients, previous.matrix_size)
    normalization = dot(reference, coefficients)
    iszero(normalization) && return nothing
    coefficients ./= normalization
    value = value_seed
    size = previous.matrix_size
    factorization = nothing
    jacobian_norm = 0.0

    for _ in 1:max_iterations
        residual_vector = matrix * coefficients - value * coefficients
        normalization_error = dot(reference, coefficients) - 1
        residual = max(norm(residual_vector), abs(normalization_error))
        # Converged once the residual is down to a small multiple of roundoff on the eigenvalue
        residual <= 100eps(Float64) * max(1.0, abs(value)) && break

        jacobian = _bordered_jacobian(matrix, value, coefficients, reference)
        factorization = lu(jacobian; check=false)
        issuccess(factorization) || return nothing
        jacobian_norm = opnorm(jacobian, 1)
        correction = factorization \ (-vcat(residual_vector, normalization_error))
        all(isfinite, correction) || return nothing
        coefficients .+= view(correction, 1:size)
        value += correction[size + 1]
    end

    #=
    The bordered Jacobian becomes singular as a neighbouring eigenvalue approaches,
    so its reciprocal condition number (cheap, given the LU factors) estimates the
    relative gap without an eigendecomposition. It underestimates the true gap, so a
    threshold on it errs toward rejecting, and the dense fallback then measures it.
    =#
    if factorization === nothing # the prediction had already converged
        jacobian = _bordered_jacobian(matrix, value, coefficients, reference)
        factorization = lu(jacobian; check=false)
        jacobian_norm = opnorm(jacobian, 1)
    end
    gap = issuccess(factorization) ? LinearAlgebra.LAPACK.gecon!(
        '1', copy(factorization.factors), jacobian_norm) : 0.0

    coefficients ./= norm(coefficients)
    phase, overlap = _phase_align!(coefficients, reference)
    angular_sep = ComplexF64(value)
    lambda = angular_sep + target^2 - 2previous.m * target
    residual = _relative_eigen_residual(
        matrix, angular_sep, coefficients)
    condition = _complex_symmetric_eigenvector_condition(coefficients)
    finite = isfinite(lambda) &&
        isfinite(residual) && isfinite(overlap) && isfinite(condition)
    finite || return nothing
    # Same :near_collision rule as in _angular_candidate
    status = gap <= sqrt(eps(Float64)) ||
        condition >= sqrt(SWSH_DEFAULT_EP_CONDITION_LIMIT) ?
        :near_collision : :predictor_corrected
    # Newton does not compare against the other eigenvectors, so there is no
    # overlap margin to report; the gap estimate guards against a close neighbour
    return AngularEigenpair(
        previous.s, previous.l, previous.m, target, previous.sheet_id,
        angular_sep, lambda, coefficients, previous.matrix_size,
        previous.truncation_order, 53, :parallel_transport,
        residual, overlap, 1.0, gap, condition, phase, status)
end

@inline function _angular_candidate_accepted(
    candidate::AngularEigenpair{Float64};
    overlap_min::Float64,
    overlap_margin_min::Float64,
    residual_tol::Float64,
)
    return all(isfinite, (candidate.lambda, candidate.residual, candidate.previous_overlap)) &&
        candidate.residual <= residual_tol &&
        candidate.previous_overlap >= overlap_min &&
        candidate.eigenvector_condition <= SWSH_DEFAULT_EP_CONDITION_LIMIT &&
        candidate.spectral_gap >= _eigenvector_gap_tol(Float64) &&
        # The lead over the runner-up matters only when the gap is small enough for roundoff to mix
        # the two eigenvectors (the :near_collision threshold); a Newton candidate reports a margin
        # of 1, so the test only bites for a dense candidate
        (candidate.spectral_gap > sqrt(eps(Float64)) ||
            candidate.overlap_margin >= overlap_margin_min)
end

function _advance_angular_segment!(
    states::Vector{AngularEigenpair{Float64}},
    previous::AngularEigenpair{Float64},
    target::ComplexF64;
    overlap_min::Float64,
    overlap_margin_min::Float64,
    residual_tol::Float64,
    min_step::Float64,
    backend::Symbol,
    depth::Int=0,
    max_depth::Int=64,
    stop_unresolved::Bool=false,
)
    target == previous.c && return previous
    if backend != :dense
        value_seed, coefficient_seed = _angular_predictor(
            states, previous, target)
        corrected = _corrected_angular_candidate(
            previous, target, value_seed, coefficient_seed)
        if corrected !== nothing && _angular_candidate_accepted(corrected;
                overlap_min, overlap_margin_min, residual_tol)
            push!(states, corrected)
            return corrected
        end
    end
    candidate = _angular_candidate(previous, target)
    if _angular_candidate_accepted(candidate;
            overlap_min, overlap_margin_min, residual_tol)
        push!(states, candidate)
        return candidate
    end

    # candidate comes from a full eigendecomposition, so its gap is exact.
    # A smaller step cannot help when the mode is degenerate at the target itself.
    if candidate.spectral_gap < _eigenvector_gap_tol(Float64)
        # A wider working precision takes the path over from here (see _advance_in_precision!)
        stop_unresolved && return nothing
        throw(AngularContinuationError(
            _unresolvable_eigenvector_message(
                previous.s, previous.l, previous.m, target,
                candidate.spectral_gap),
            previous.c, target, candidate.previous_overlap,
            candidate.overlap_margin, candidate.residual,
            candidate.spectral_gap, candidate.eigenvector_condition))
    end

    step = abs(target - previous.c)
    scale = max(1.0, abs(previous.c), abs(target))
    if depth >= max_depth || step <= min_step * scale
        throw(AngularContinuationError(
            "SWSH angular continuation exhausted its minimum step (or its $max_depth halvings); " *
            "the path may cross an eigenvalue collision or exceptional point.",
            previous.c, target, candidate.previous_overlap,
            candidate.overlap_margin, candidate.residual,
            candidate.spectral_gap, candidate.eigenvector_condition))
    end

    midpoint = previous.c + (target - previous.c) / 2
    middle = _advance_angular_segment!(
        states, previous, midpoint;
        overlap_min, overlap_margin_min, residual_tol, min_step, backend,
        depth=depth + 1, max_depth, stop_unresolved)
    middle === nothing && return nothing
    return _advance_angular_segment!(
        states, middle, target;
        overlap_min, overlap_margin_min, residual_tol, min_step, backend,
        depth=depth + 1, max_depth, stop_unresolved)
end

#=
One step of the path in the working precision T, where Float64 can no longer tell the mode from a
neighbour: Rayleigh quotient iteration on the band λ-matrix at the target, started from the
extrapolation of the previous states. The step is accepted when the result stays on the mode (its
overlap with the previous state) and has converged; otherwise it is halved, down to `min_step`.
The gap to the neighbours is not measured here (that needs an eigendecomposition in T), so the
states carry NaN for it and for the overlap margin.
=#
function _advance_in_precision!(
    states::Vector{AngularEigenpair{T}},
    previous::AngularEigenpair{T},
    target::Complex{T};
    overlap_min::Real,
    residual_tol::Real,
    min_step::Real,
    depth::Int=0,
    max_depth::Int=64,
) where {T}
    target == previous.c && return previous
    value_seed, coefficient_seed = _angular_predictor(states, previous, target)
    band = _spectral_band(target, previous.s, previous.m, previous.matrix_size)
    shift = target^2 - 2previous.m * target
    lambda, coefficients, residual = _rayleigh_quotient_iteration(
        band, value_seed + shift, coefficient_seed)
    coefficients = coefficients / norm(coefficients)
    phase, overlap = _phase_align!(coefficients, previous.coefficients)
    condition = _complex_symmetric_eigenvector_condition(coefficients)
    if overlap >= overlap_min && residual <= _tol(T, residual_tol)
        state = AngularEigenpair(
            previous.s, previous.l, previous.m, target, previous.sheet_id,
            lambda - shift, lambda, coefficients, previous.matrix_size,
            previous.truncation_order, precision(T), :parallel_transport,
            residual, overlap, T(NaN), T(NaN), condition, phase, :precision_tracked)
        push!(states, state)
        return state
    end
    step = abs(target - previous.c)
    scale = max(one(T), abs(previous.c), abs(target))
    if depth >= max_depth || step <= min_step * scale
        throw(AngularContinuationError(
            "SWSH angular continuation in $(_precision_name(T)) exhausted its minimum step " *
            "(or its $max_depth halvings); " *
            "the path may cross an eigenvalue collision or exceptional point.",
            ComplexF64(previous.c), ComplexF64(target), Float64(overlap), NaN,
            Float64(residual), NaN, Float64(condition)))
    end
    midpoint = previous.c + (target - previous.c) / 2
    middle = _advance_in_precision!(states, previous, midpoint;
        overlap_min, residual_tol, min_step, depth=depth + 1, max_depth)
    return _advance_in_precision!(states, middle, target;
        overlap_min, residual_tol, min_step, depth=depth + 1, max_depth)
end

#=
Along a real path the banded solver needs no tracking: for real c the eigenvalues
keep their order, so the mode is simply the one at its index. Each point is solved
on its own, then phase-aligned with the previous point so that the states follow
the same convention as the other backends.
=#
function _banded_angular_state(
    previous::AngularEigenpair{Float64},
    target::ComplexF64,
)
    c = real(target)
    # Throws if the eigenvector cannot be resolved (see _eigenvector_gap_tol)
    lambda, coefficients, gap = _real_lambda_eigenpair_at_size(
        c, previous.s, previous.l, previous.m, previous.matrix_size)
    phase, overlap = _phase_align!(coefficients, previous.coefficients)
    angular_sep = lambda - muladd(c, c, -2previous.m * c)
    residual = _real_lambda_residual(
        c, previous.s, previous.m, lambda, coefficients)
    condition = _complex_symmetric_eigenvector_condition(coefficients)
    # The mode is picked by index, not by overlap, so there is no overlap margin
    return AngularEigenpair(
        previous.s, previous.l, previous.m, target, previous.sheet_id,
        ComplexF64(angular_sep), ComplexF64(lambda), coefficients,
        previous.matrix_size, previous.truncation_order, 53,
        :parallel_transport, residual, overlap, 1.0, gap, condition,
        ComplexF64(phase), :banded)
end

@doc raw"""
    track_angular_mode(s::Int, l::Int, m::Int, requested_path; backend=:auto, sheet_id::Symbol=:principal, truncation_order=nothing, max_step::Real=0.5, min_step::Real=1e-10, overlap_min::Real=0.65, overlap_margin_min::Real=1e-8, residual_tol::Real=5e-12, precision=nothing)

Follow the spin-weighted spheroidal mode of spin weight `s`, harmonic index `l`, and
azimuthal index `m` along `requested_path`, a sequence of spheroidicities `c` ($c = a\omega$)
visited in order. The path always starts at `c = 0`, where the mode is exactly the
spin-weighted spherical harmonic of index `l`; `0` is prepended if the path does not
start there. Each leg of the path is divided into steps of at most `max_step`.

Because the mode is followed continuously from `c = 0`, "mode `l`" means the mode
connected to the spherical harmonic `l` along *this* path. For complex `c` the result
can depend on the path once it passes a branch point of the eigenvalue: a different
path can end on a different mode. `spin_weighted_spheroidal_harmonic` and
`spin_weighted_spheroidal_eigenvalue` use the straight line from `0` to `c`.

The `backend` argument controls how each step is solved:
- `:auto` (default): Newton iteration seeded by extrapolating the previous steps, falling
  back to a full dense eigendecomposition when Newton fails or its result is rejected,
- `:dense`: a full dense eigendecomposition at every step, keeping the eigenvector with
  the largest overlap with the previous step,
- `:banded`: real paths only. Every point is solved on its own with the banded eigenpair
  solver, which picks the mode by its position among the (real) eigenvalues; for real `c`
  this is the same mode as following it continuously. The matrix size is fixed by
  `truncation_order` (it is not increased adaptively), and the acceptance tolerances
  below are not used. A path with any complex point throws an `ArgumentError`.
A backend may also be given as a string; names are case-insensitive.

The other keyword arguments are:
- `sheet_id`: a label stored in every returned state, to tell paths apart,
- `truncation_order`: the number of spherical harmonics kept beyond the target mode,
  so that the matrix has `l - max(|m|, |s|) + 1 + truncation_order` rows. The size is
  fixed along the path; by default (`nothing`) it is chosen for the largest `|c|` on the
  path, growing like `sqrt(|c|)` so that the neglected coefficients stay below about `1e-14`,
- `min_step`: when a step is rejected it is halved; continuation gives up once a step
  is smaller than `min_step * max(1, |c|)`,
- `overlap_min`: the smallest overlap with the previous step's eigenvector for a step to be accepted,
- `overlap_margin_min`: when the mode is close to another eigenvalue, how much larger its
  overlap must be than the runner-up's,
- `residual_tol`: the largest relative residual of the eigen-equation for a step to be accepted,
- `precision`: compute in BigFloat at this many bits (the path is converted once); by default the
  path is followed in its own floating-point type (at least Float64).

The path is followed in Float64, and in a wider precision every state is then refined at its own
`c`, so the acceptance tolerances above apply to the Float64 path.

Return an `AngularPathResult` whose `states` hold an `AngularEigenpair` for `c = 0`, every
intermediate step, and every point of the path, in order; `last(result.states)` is the mode
at the end of the path, and can be passed to `spin_weighted_spheroidal_harmonic`. When the
path ends where it started, `monodromy_overlap` and `lambda_closure_error` compare the final
state with the initial one, and `status` is `:closed` or `:monodromy_failed`; otherwise
`status` is `:open`.

Throw an `AngularContinuationError` when a step cannot be accepted even at the smallest
step size (or after 64 halvings), or when the mode is degenerate with a neighbouring one at
the end point of a step, beyond what the working precision resolves; either usually means
the path runs into an eigenvalue collision or an exceptional point. With `backend=:banded`
the degenerate case is reported by the banded solver as an `ErrorException` instead. An
empty path throws an `ArgumentError`.

# Example
```julia
path = [0.0, 0.5 - 2.0im, 1.0 - 4.0im]
result = track_angular_mode(-2, 2, 2, path)
pair = last(result.states)
S = spin_weighted_spheroidal_harmonic(pair)
S(1.1, 0.3)
```
"""
function track_angular_mode(
    s::Int,
    l::Int,
    m::Int,
    requested_path;
    precision::Union{Nothing,Integer}=nothing,
    kwargs...,
)
    _check_precision_keyword(precision)
    path = _with_precision_keyword(collect(requested_path), precision)
    isempty(path) && throw(ArgumentError("requested_path must not be empty."))
    T = _float_type(path...)
    bits = precision === nothing ? _input_precision(path) : precision
    return _with_precision(T, bits) do
        _track_angular_mode(T, s, l, m, Complex{T}.(path); kwargs...)
    end
end

#=
The path is followed in Float64; in a wider working precision the states are then refined one by
one at their own c (all of them, or only the last one when that is all the caller needs), with the
points of the requested path taken exactly in the working precision.
=#
function _track_angular_mode(
    ::Type{T},
    s::Int,
    l::Int,
    m::Int,
    path::Vector{Complex{T}};
    backend=:auto,
    sheet_id::Symbol=:principal,
    truncation_order::Union{Int,Nothing}=nothing,
    max_step::Real=SWSH_DEFAULT_MAX_PATH_STEP,
    min_step::Real=SWSH_DEFAULT_MIN_PATH_STEP,
    overlap_min::Real=SWSH_DEFAULT_OVERLAP_MIN,
    overlap_margin_min::Real=SWSH_DEFAULT_OVERLAP_MARGIN_MIN,
    residual_tol::Real=SWSH_DEFAULT_RESIDUAL_TOL,
    refine_all::Bool=true,
) where {T}
    selected_backend = _spectral_backend(backend)
    if selected_backend == :banded && !all(isreal, path)
        throw(ArgumentError(
            "backend=:banded requires real c, so every point on the path " *
            "must be real; use :auto or :dense."))
    end
    first(path) == 0 || pushfirst!(path, zero(Complex{T}))
    if truncation_order === nothing
        # The matrix size is fixed along the path, so size it for the farthest point
        truncation_order = _complex_truncation_order(
            s, l, m, -1, maximum(abs, path))
    end
    matrix_size = _angular_matrix_size(s, l, m, truncation_order)
    current = _spherical_angular_pair(
        s, l, m, matrix_size, truncation_order, sheet_id)
    states = AngularEigenpair{Float64}[current]
    anchors = [1] # the state at each point of the path
    # Where the Float64 path stops because it cannot resolve the mode (only in a wider precision,
    # which takes it over from there): the index of the requested point it was heading for
    resume = 0

    for (k, requested_target) in enumerate(ComplexF64.(path))
        k == 1 && continue
        segment = requested_target - current.c
        subdivisions = max(1, ceil(Int, abs(segment) / Float64(max_step)))
        for step_index in 1:subdivisions
            target = current.c +
                (requested_target - current.c) / (subdivisions - step_index + 1)
            if selected_backend == :banded
                # Same steps as the other backends, so their states line up
                current = _banded_angular_state(current, target)
                push!(states, current)
            else
                advanced = _advance_angular_segment!(
                    states, current, target;
                    overlap_min=Float64(overlap_min),
                    overlap_margin_min=Float64(overlap_margin_min),
                    residual_tol=Float64(residual_tol),
                    min_step=Float64(min_step), backend=selected_backend,
                    stop_unresolved=!_seed_is_final(T))
                advanced === nothing && (resume = k; break)
                current = advanced
            end
        end
        resume > 0 && break
        push!(anchors, length(states))
    end

    states = if _seed_is_final(T)
        states
    else
        # Every state is refined, each phase-aligned with the refined state before it, so that the
        # parallel transport along the path is carried out in T
        exact = Complex{T}.(getfield.(states, :c))
        exact[anchors] .= path[1:length(anchors)]
        refined = AngularEigenpair{T}[_refine_angular_pair(first(states), first(exact), nothing)]
        for k in 2:length(states)
            push!(refined, _refine_angular_pair(states[k], exact[k], last(refined).coefficients))
        end
        refine_all || resume > 0 ? refined : refined[end:end]
    end
    if resume > 0
        for requested_target in path[resume:end]
            previous = last(states)
            subdivisions = max(1, ceil(Int, abs(requested_target - previous.c) / max_step))
            for step_index in 1:subdivisions
                target = previous.c + (requested_target - previous.c) / (subdivisions - step_index + 1)
                previous = _advance_in_precision!(states, previous, target;
                    overlap_min, residual_tol, min_step)
            end
        end
        refine_all || (states = states[end:end])
    end

    # The path is closed when its last point returns to its first (always c = 0) up to a few units
    # of roundoff
    closed = length(path) > 1 &&
        abs(last(path) - first(path)) <=
        16eps(T) * max(one(T), abs(first(path)))
    if closed && refine_all
        initial = first(states)
        final = last(states)
        initial_coefficients = _resize_coefficients(
            initial.coefficients, length(final.coefficients))
        monodromy = dot(initial_coefficients, final.coefficients) /
            max(norm(initial_coefficients) * norm(final.coefficients),
                floatmin(T))
        closure_error = abs(final.lambda - initial.lambda) /
            max(one(T), abs(initial.lambda), abs(final.lambda))
        # On the same mode the final eigenvalue agrees with the initial (spherical) one to roundoff; on
        # another mode, which happens when the loop encircles a branch point and the modes are
        # exchanged, the two differ by far more, so the threshold only has to separate the two cases
        status = closure_error <= _tol(T, 1.0e-10) ? :closed : :monodromy_failed
    else
        monodromy = Complex{T}(NaN, NaN)
        closure_error = T(NaN)
        status = :open
    end
    return AngularPathResult(
        states, path, sheet_id, monodromy, closure_error, status)
end

@doc raw"""
    continue_angular_mode(s::Int, l::Int, m::Int, c; backend=:auto, sheet_id::Symbol=:straight_from_spherical, truncation_order=nothing, cache=DEFAULT_ANGULAR_CACHE, precision=nothing, kwargs...)

Follow the mode `l` along the straight path from `c = 0` to `c` (see `track_angular_mode`, which
takes the remaining keyword arguments) and return its `AngularEigenpair` at `c`. The result is
computed in the floating-point type of `c` (at least Float64), or in BigFloat at `precision` bits.
Results are kept in `cache` (pass `nothing` to bypass it), keyed by everything they depend on:
the mode, `c`, its type and precision, the truncation, the backend and the keyword arguments.
"""
function continue_angular_mode(
    s::Int,
    l::Int,
    m::Int,
    c;
    backend=:auto,
    sheet_id::Symbol=:straight_from_spherical,
    truncation_order::Union{Int,Nothing}=nothing,
    cache::Union{AngularCache,Nothing}=DEFAULT_ANGULAR_CACHE,
    precision::Union{Nothing,Integer}=nothing,
    kwargs...,
)
    _check_precision_keyword(precision)
    c = _with_precision_keyword(c, precision)
    T = _float_type(c)
    bits = precision === nothing ? _input_precision(c) : precision
    return _with_precision(T, bits) do
        target = Complex{T}(c)
        selected_backend = _spectral_backend(backend)
        # Resolved before the cache lookup, so that the key holds the size actually used
        truncation_order === nothing &&
            (truncation_order = _complex_truncation_order(s, l, m, -1, target))
        key = AngularCacheKey(s, l, m, target, sheet_id, Base.precision(T),
            truncation_order, selected_backend, NamedTuple(kwargs))
        if cache !== nothing
            cached = _angular_cache_get(cache, key)
            cached === nothing || return cached
        end
        # At most 64 requested points: track_angular_mode subdivides each leg further to max_step,
        # so this bounds only the length of the requested path, not the step
        steps = max(1, min(64,
            ceil(Int, abs(ComplexF64(target)) / SWSH_DEFAULT_MAX_PATH_STEP)))
        path = Complex{T}[
            target * (index / steps) for index in 0:steps
        ]
        result = _track_angular_mode(T,
            s, l, m, path;
            sheet_id, truncation_order, backend=selected_backend,
            refine_all=false, kwargs...)
        pair = last(result.states)
        cache === nothing ? pair : _angular_cache_put!(cache, key, pair)
    end
end

"""
    continue_angular_lateral_pair(s, l, m, a, sigma, epsilon; truncation_order=nothing, precision=nothing, kwargs...)

Follow the mode `l` to the two points `c = ±a epsilon - i a sigma` on either side of the imaginary axis,
along the paths `0 → -i a sigma → c` (with `sigma > 0`, `epsilon > 0`). Return a named tuple `(right, left)`
of the two `AngularPathResult`s. The remaining keyword arguments are those of `track_angular_mode`.
"""
function continue_angular_lateral_pair(
    s::Int,
    l::Int,
    m::Int,
    a::Real,
    sigma::Real,
    epsilon::Real;
    truncation_order::Union{Int,Nothing}=nothing,
    precision::Union{Nothing,Integer}=nothing,
    kwargs...,
)
    _check_precision_keyword(precision)
    sigma > 0 || throw(ArgumentError("sigma must be positive."))
    epsilon > 0 || throw(ArgumentError("epsilon must be positive."))
    a, sigma, epsilon = (_with_precision_keyword(x, precision) for x in (a, sigma, epsilon))
    T = _float_type(a, sigma, epsilon)
    bits = precision === nothing ? _input_precision(a, sigma, epsilon) : precision
    return _with_precision(T, bits) do
        common = Complex{T}(zero(T), -T(a) * T(sigma))
        right_target = Complex{T}(
            T(a) * T(epsilon), imag(common))
        left_target = Complex{T}(
            -T(a) * T(epsilon), imag(common))
        right = _track_angular_mode(T,
            s, l, m, Complex{T}[0, common, right_target];
            sheet_id=:right_lateral, truncation_order, kwargs...)
        left = _track_angular_mode(T,
            s, l, m, Complex{T}[0, common, left_target];
            sheet_id=:left_lateral, truncation_order, kwargs...)
        (right=right, left=left)
    end
end

"""
    mirror_angular_parameters(s, l, m, c)

The parameters `(s, l, -m, -conj(c))` of the mirror mode, whose eigenvalue is the complex conjugate of
that of `(s, l, m, c)`.
"""
@inline function mirror_angular_parameters(
    s::Int,
    l::Int,
    m::Int,
    c,
)
    return (s=s, l=l, m=-m, c=-conj(c))
end

"""
    angular_mirror_residual(pair, mirrored)

The relative difference between the eigenvalue of `mirrored`, the eigenpair at
`mirror_angular_parameters` of `pair`, and the complex conjugate of the eigenvalue of `pair`. It is a
check on both: it vanishes up to roundoff.
"""
function angular_mirror_residual(
    pair::AngularEigenpair,
    mirrored::AngularEigenpair,
)
    expected = mirror_angular_parameters(
        pair.s, pair.l, pair.m, pair.c)
    mirrored.s == expected.s && mirrored.l == expected.l &&
        mirrored.m == expected.m ||
        throw(ArgumentError("mirrored pair has inconsistent mode indices."))
    isapprox(mirrored.c, expected.c; rtol=0, atol=0) ||
        throw(ArgumentError("mirrored pair has inconsistent spheroidicity."))
    return abs(mirrored.lambda - conj(pair.lambda)) /
        max(one(real(typeof(pair.lambda))), abs(mirrored.lambda), abs(pair.lambda))
end

#=
Refine a Float64 state at c in the working precision of c: Rayleigh quotient iteration on the band
λ-matrix from the Float64 eigenpair, phase-aligned with the refined previous state of the path
(`previous`, the coefficients; `nothing` for the first state, whose phase is that of the seed). The
diagnostics of the path (overlaps, gap) are those of the Float64 state.
=#
function _refine_angular_pair(seed::AngularEigenpair{Float64}, c::Complex{T}, previous) where {T}
    band = _spectral_band(c, seed.s, seed.m, seed.matrix_size)
    lambda, coefficients, residual = _rayleigh_quotient_iteration(
        band, Complex{T}(seed.lambda), Complex{T}.(seed.coefficients))
    coefficients = coefficients / norm(coefficients)
    phase, _ = _phase_align!(coefficients, something(previous, Complex{T}.(seed.coefficients)))
    angular_sep = lambda - (c^2 - 2seed.m * c)
    condition = _complex_symmetric_eigenvector_condition(coefficients)
    return AngularEigenpair(
        seed.s, seed.l, seed.m, c, seed.sheet_id, angular_sep, lambda,
        coefficients, seed.matrix_size, seed.truncation_order,
        precision(T), seed.normalization, residual, T(seed.previous_overlap),
        T(seed.overlap_margin), T(seed.spectral_gap), condition, Complex{T}(seed.phase_factor) * phase,
        seed.status)
end

"""
    angular_precision_certificate(s, l, m, c; truncation_orders=(24, 32), precision_bits=(128, 160), kwargs...)

Compute the mode `l` at `c` with `continue_angular_mode` at two truncation orders, and at the higher one
also at two `BigFloat` precisions, and compare them. Return a named tuple with `accepted` (the eigenvalue
is stable to `1e-12` under truncation and to `1e-14` under precision, the residual is small and the two
eigenvectors agree), the four eigenpairs and the drifts, overlap, residual and eigenvector condition.
"""
function angular_precision_certificate(
    s::Int,
    l::Int,
    m::Int,
    c;
    sheet_id::Symbol=:straight_from_spherical,
    truncation_orders::Tuple{Int,Int}=(24, 32),
    precision_bits::Tuple{Int,Int}=(128, 160),
    kwargs...,
)
    low_order, high_order = truncation_orders
    low_bits, high_bits = precision_bits
    low_pair = continue_angular_mode(
        s, l, m, c;
        sheet_id=Symbol(sheet_id, :_order_low),
        truncation_order=low_order, cache=nothing, kwargs...)
    high_pair, refined_low, refined_high = (continue_angular_mode(
        s, l, m, c;
        sheet_id=Symbol(sheet_id, :_order_high),
        truncation_order=high_order, cache=nothing, precision=bits, kwargs...)
        for bits in (nothing, low_bits, high_bits))
    low_coefficients = _resize_coefficients(
        low_pair.coefficients, length(high_pair.coefficients))
    R = real(eltype(high_pair.coefficients))
    vector_overlap = abs(dot(low_coefficients, high_pair.coefficients)) /
        max(norm(low_coefficients) * norm(high_pair.coefficients),
            floatmin(R))
    truncation_drift = abs(low_pair.lambda - high_pair.lambda) /
        max(one(R), abs(low_pair.lambda), abs(high_pair.lambda))
    precision_drift = Float64(abs(refined_low.lambda - refined_high.lambda) /
        max(one(BigFloat), abs(refined_low.lambda), abs(refined_high.lambda)))
    # Accept when the truncation error (the drift between the two orders) and the Float64 residual are
    # both at the 1e-12 level, the eigenvectors of the two orders agree to 1e-10, and the two
    # refinements agree far more closely: they converge to the same eigenvalue of the same matrix, so a
    # precision drift above 1e-14 means a refinement did not converge
    accepted = truncation_drift <= 1.0e-12 &&
        precision_drift <= 1.0e-14 &&
        high_pair.residual <= SWSH_DEFAULT_RESIDUAL_TOL &&
        vector_overlap >= 1 - 1.0e-10
    return (
        accepted=accepted,
        low_order=low_pair,
        high_order=high_pair,
        low_precision=refined_low,
        high_precision=refined_high,
        truncation_drift=truncation_drift,
        precision_drift=precision_drift,
        eigenvector_overlap=vector_overlap,
        residual=high_pair.residual,
        eigenvector_condition=high_pair.eigenvector_condition,
    )
end


#=
The overall phase of a harmonic: its overlap with sYlm, the coefficient of the target mode, real and
positive, as for real c and in Leaver's method. The states of a path keep the phase of the parallel
transport, which is what keeps the path continuous; a harmonic takes its phase from the target
coefficient instead. When that overlap is too small to define a phase (below 1e-3 of the norm, the
threshold of Leaver's method) the phase is left as it is.
=#
function _overlap_phase_fixed(coefficients::AbstractVector, s::Int, l::Int, m::Int)
    b = coefficients[l - max(abs(s), abs(m)) + 1]
    (isreal(b) && real(b) > 0) && return coefficients
    abs(b) > 1e-3 * norm(coefficients) || return coefficients
    return coefficients .* (abs(b) / b)
end

"""
    spin_weighted_spheroidal_harmonic(pair::AngularEigenpair; method="auto", normalization=:modulus)

Build the harmonic from a continued eigenpair (see `continue_angular_mode` and `track_angular_mode`),
in the floating-point type of the pair. `method` is `"auto"`, `"direct"`, `"jacobi"` or `"chebyshev"`.
Like every harmonic, it is normalized with its modulus and its phase is fixed by making the overlap
with the spherical harmonic `l` real and positive, whatever the phase of the pair. With
`normalization=:analytic` it is instead normalized so that ∫S² dΩ = 1, the analytic continuation of
real `c`; see `spin_weighted_spheroidal_harmonic(s, l, m, c)`.
"""
function spin_weighted_spheroidal_harmonic(
    pair::AngularEigenpair{T};
    method="auto",
    normalization=:modulus,
) where {T}
    formatted_method = _format_method_name(method)
    analytic = _format_normalization(normalization) == :analytic
    # Leaver's method finds its own eigenpair, so it has no use for this one
    formatted_method == "leaver" && throw(ArgumentError(
        "A continued eigenpair is evaluated from its spherical-harmonic coefficients, so " *
        "method must be \"auto\", \"direct\", \"jacobi\" or \"chebyshev\" here. For Leaver's " *
        "method, call spin_weighted_spheroidal_harmonic(s, l, m, c; method=\"leaver\")."))
    bits = _input_precision(pair.c)
    return _with_precision(T, bits) do
        coefficients_params = SpectralDecompositionInputParams(
            pair.s, pair.l, pair.m, pair.c, pair.matrix_size)
        l_list = construct_all_l_in_matrix(
            pair.s, pair.m, pair.matrix_size)
        normalization = one(real(eltype(pair.coefficients)))
        coefficients = _overlap_phase_fixed(pair.coefficients, pair.s, pair.l, pair.m)
        if analytic
            coefficients = _analytically_normalized(coefficients, pair.s, pair.l, pair.m)
        end
        if formatted_method == "chebyshev"
            # The eigenpair supplies lambda, and through its coefficients the values at the equator
            S_half, dS_half = _spheroidal_equator_values(T,
                coefficients_params, coefficients)
            angular_sep = pair.lambda - pair.c^2 + 2 * pair.m * pair.c
            chebyshev_solution = _solve_spheroidal_harmonic_chebyshev(T,
                pair.s, pair.l, pair.m, pair.c, angular_sep, S_half, dS_half)
            spherical_harmonics_l = Vector{Union{
                SpinWeightedSphericalHarmonicFunction, Nothing}}(
                nothing, length(l_list))
            return SpinWeightedSpheroidalHarmonicFunction{T}(
                coefficients_params, coefficients, spherical_harmonics_l,
                normalization, pair.lambda, :chebyshev, chebyshev_solution, nothing, bits)
        end
        spherical_harmonics_l = [
            abs(coefficients[index]) >= _spherical_coefficient_tol(T) ?
            spin_weighted_spherical_harmonic(
                pair.s, l_list[index], pair.m; method=formatted_method,
                precision=_precision_keyword(bits)) :
            nothing for index in eachindex(l_list)
        ]
        SpinWeightedSpheroidalHarmonicFunction{T}(
            coefficients_params, coefficients, spherical_harmonics_l,
            normalization, pair.lambda, :spectral, nothing, nothing, bits)
    end
end

"""
    angular_observables(pair::AngularEigenpair, theta, phi=0; derivative_order=2, method="auto")

The value of the harmonic of `pair` at `(theta, phi)` and its `theta`-derivatives up to
`derivative_order` (`0`, `1` or `2`), with `lambda`, `sheet_id` and `residual` of the pair, as a named
tuple `(value, first_derivative, second_derivative, lambda, sheet_id, residual)`.
"""
function angular_observables(
    pair::AngularEigenpair,
    theta,
    phi=0;
    derivative_order::Int=2,
    method="auto",
)
    0 <= derivative_order <= 2 ||
        throw(ArgumentError("derivative_order must be 0, 1, or 2."))
    harmonic = spin_weighted_spheroidal_harmonic(pair; method)
    value = harmonic(theta, phi)
    first_derivative = derivative_order >= 1 ?
        harmonic(theta, phi; theta_derivative=1) : nothing
    second_derivative = derivative_order >= 2 ?
        harmonic(theta, phi; theta_derivative=2) : nothing
    return (
        value=value,
        first_derivative=first_derivative,
        second_derivative=second_derivative,
        lambda=pair.lambda,
        sheet_id=pair.sheet_id,
        residual=pair.residual,
    )
end
