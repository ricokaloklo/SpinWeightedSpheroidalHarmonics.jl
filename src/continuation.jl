const SWSH_DEFAULT_ANGULAR_ORDER = 32
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
end

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

struct AngularPathResult{T<:AbstractFloat}
    states::Vector{AngularEigenpair{T}}
    requested_path::Vector{Complex{T}}
    sheet_id::Symbol
    monodromy_overlap::Complex{T}
    lambda_closure_error::T
    status::Symbol
end

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

@inline _binary64_angular_input(::Float16) = true
@inline _binary64_angular_input(::Float32) = true
@inline _binary64_angular_input(::Float64) = true
@inline _binary64_angular_input(value::Complex) =
    _binary64_angular_input(real(value)) &&
    _binary64_angular_input(imag(value))
@inline _binary64_angular_input(value) = false

function _require_binary64_angular_input(c)
    _binary64_angular_input(c) && return c
    throw(ArgumentError(
        "continued SWSH currently accepts binary64-compatible c only; " *
        "use angular_precision_certificate for 128/160-bit refinement."))
end

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

function _relative_spectral_gap(values, index::Int)
    length(values) == 1 && return Inf
    value = values[index]
    gap = minimum(abs(values[j] - value) for j in eachindex(values)
        if j != index)
    return gap / max(one(real(value)), abs(value))
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
    gap = _relative_spectral_gap(decomposition.values, selected)
    condition = _complex_symmetric_eigenvector_condition(coefficients)
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
    states::Vector{AngularEigenpair{Float64}},
    previous::AngularEigenpair{Float64},
    target::ComplexF64,
)
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

    for _ in 1:max_iterations
        residual_vector = matrix * coefficients - value * coefficients
        normalization_error = dot(reference, coefficients) - 1
        residual = max(norm(residual_vector), abs(normalization_error))
        residual <= 100eps(Float64) * max(1.0, abs(value)) && break

        jacobian = zeros(ComplexF64, size + 1, size + 1)
        jacobian[1:size, 1:size] .= matrix
        @inbounds for index in 1:size
            jacobian[index, index] -= value
            jacobian[index, size + 1] = -coefficients[index]
            jacobian[size + 1, index] = conj(reference[index])
        end
        factorization = lu(jacobian; check=false)
        issuccess(factorization) || return nothing
        correction = factorization \ (-vcat(residual_vector, normalization_error))
        all(isfinite, correction) || return nothing
        coefficients .+= view(correction, 1:size)
        value += correction[size + 1]
    end

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
    status = condition >= sqrt(SWSH_DEFAULT_EP_CONDITION_LIMIT) ?
        :near_collision : :predictor_corrected
    return AngularEigenpair(
        previous.s, previous.l, previous.m, target, previous.sheet_id,
        angular_sep, lambda, coefficients, previous.matrix_size,
        previous.truncation_order, 53, :parallel_transport,
        residual, overlap, 1.0, Inf, condition, phase, status)
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

    step = abs(target - previous.c)
    scale = max(1.0, abs(previous.c), abs(target))
    if depth >= max_depth || step <= min_step * scale
        throw(AngularContinuationError(
            "SWSH angular continuation exhausted its minimum step; " *
            "the path may cross an eigenvalue collision or exceptional point.",
            previous.c, target, candidate.previous_overlap,
            candidate.overlap_margin, candidate.residual,
            candidate.spectral_gap, candidate.eigenvector_condition))
    end

    midpoint = previous.c + (target - previous.c) / 2
    middle = _advance_angular_segment!(
        states, previous, midpoint;
        overlap_min, overlap_margin_min, residual_tol, min_step, backend,
        depth=depth + 1, max_depth)
    return _advance_angular_segment!(
        states, middle, target;
        overlap_min, overlap_margin_min, residual_tol, min_step, backend,
        depth=depth + 1, max_depth)
end

function track_angular_mode(
    s::Int,
    l::Int,
    m::Int,
    requested_path;
    backend=:auto,
    sheet_id::Symbol=:principal,
    truncation_order::Int=SWSH_DEFAULT_ANGULAR_ORDER,
    max_step::Real=SWSH_DEFAULT_MAX_PATH_STEP,
    min_step::Real=SWSH_DEFAULT_MIN_PATH_STEP,
    overlap_min::Real=SWSH_DEFAULT_OVERLAP_MIN,
    overlap_margin_min::Real=SWSH_DEFAULT_OVERLAP_MARGIN_MIN,
    residual_tol::Real=SWSH_DEFAULT_RESIDUAL_TOL,
)
    selected_backend = _spectral_backend(backend)
    path = ComplexF64.(collect(requested_path))
    isempty(path) && throw(ArgumentError("requested_path must not be empty."))
    first(path) == 0 || pushfirst!(path, 0.0 + 0.0im)
    matrix_size = _angular_matrix_size(s, l, m, truncation_order)
    current = _spherical_angular_pair(
        s, l, m, matrix_size, truncation_order, sheet_id)
    states = AngularEigenpair{Float64}[current]

    for requested_target in Iterators.drop(path, 1)
        segment = requested_target - current.c
        subdivisions = max(1, ceil(Int, abs(segment) / Float64(max_step)))
        for step_index in 1:subdivisions
            target = current.c +
                (requested_target - current.c) / (subdivisions - step_index + 1)
            current = _advance_angular_segment!(
                states, current, target;
                overlap_min=Float64(overlap_min),
                overlap_margin_min=Float64(overlap_margin_min),
                residual_tol=Float64(residual_tol),
                min_step=Float64(min_step), backend=selected_backend)
        end
    end

    closed = length(path) > 1 &&
        abs(last(path) - first(path)) <=
        16eps(Float64) * max(1.0, abs(first(path)))
    if closed
        initial = first(states)
        final = last(states)
        initial_coefficients = _resize_coefficients(
            initial.coefficients, length(final.coefficients))
        monodromy = dot(initial_coefficients, final.coefficients) /
            max(norm(initial_coefficients) * norm(final.coefficients),
                floatmin(Float64))
        closure_error = abs(final.lambda - initial.lambda) /
            max(1.0, abs(initial.lambda), abs(final.lambda))
        status = closure_error <= 1.0e-10 ? :closed : :monodromy_failed
    else
        monodromy = ComplexF64(NaN, NaN)
        closure_error = NaN
        status = :open
    end
    return AngularPathResult(
        states, path, sheet_id, monodromy, closure_error, status)
end

function continue_angular_mode(
    s::Int,
    l::Int,
    m::Int,
    c;
    backend=:auto,
    sheet_id::Symbol=:straight_from_spherical,
    truncation_order::Int=SWSH_DEFAULT_ANGULAR_ORDER,
    cache::Union{AngularCache,Nothing}=DEFAULT_ANGULAR_CACHE,
    kwargs...,
)
    _require_binary64_angular_input(c)
    selected_backend = _spectral_backend(backend)
    c64 = ComplexF64(c)
    key = AngularCacheKey(
        s, l, m, c64, sheet_id, 53, truncation_order, selected_backend)
    if cache !== nothing
        cached = _angular_cache_get(cache, key)
        cached === nothing || return cached
    end
    steps = max(1, min(64,
        ceil(Int, abs(c64) / SWSH_DEFAULT_MAX_PATH_STEP)))
    path = ComplexF64[
        c64 * (index / steps) for index in 0:steps
    ]
    result = track_angular_mode(
        s, l, m, path;
        sheet_id, truncation_order, backend=selected_backend, kwargs...)
    pair = last(result.states)
    return cache === nothing ? pair : _angular_cache_put!(cache, key, pair)
end

function continue_angular_lateral_pair(
    s::Int,
    l::Int,
    m::Int,
    a::Real,
    sigma::Real,
    epsilon::Real;
    truncation_order::Int=SWSH_DEFAULT_ANGULAR_ORDER,
    kwargs...,
)
    sigma > 0 || throw(ArgumentError("sigma must be positive."))
    epsilon > 0 || throw(ArgumentError("epsilon must be positive."))
    common = ComplexF64(0.0, -Float64(a) * Float64(sigma))
    right_target = ComplexF64(
        Float64(a) * Float64(epsilon), imag(common))
    left_target = ComplexF64(
        -Float64(a) * Float64(epsilon), imag(common))
    right = track_angular_mode(
        s, l, m, ComplexF64[0, common, right_target];
        sheet_id=:right_lateral, truncation_order, kwargs...)
    left = track_angular_mode(
        s, l, m, ComplexF64[0, common, left_target];
        sheet_id=:left_lateral, truncation_order, kwargs...)
    return (right=right, left=left)
end

@inline function mirror_angular_parameters(
    s::Int,
    l::Int,
    m::Int,
    c,
)
    return (s=s, l=l, m=-m, c=-conj(c))
end

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
        max(1.0, abs(mirrored.lambda), abs(pair.lambda))
end

function _big_swsh_f(::Type{T}, s::Int, l::Int, m::Int) where {T<:AbstractFloat}
    (l == -1 && iszero(m) && iszero(s)) && return zero(T)
    lp1 = T(l + 1)
    return sqrt(((lp1^2 - T(m)^2) /
        (T(2l + 3) * T(2l + 1))) *
        ((lp1^2 - T(s)^2) / lp1^2))
end

function _big_swsh_g(::Type{T}, s::Int, l::Int, m::Int) where {T<:AbstractFloat}
    iszero(l) && return zero(T)
    ll = T(l)
    return sqrt(((ll^2 - T(m)^2) / (T(4) * ll^2 - one(T))) *
        ((ll^2 - T(s)^2) / ll^2))
end

function _big_swsh_h(::Type{T}, s::Int, l::Int, m::Int) where {T<:AbstractFloat}
    (iszero(l) || iszero(s)) && return zero(T)
    return -T(m * s) / T(l * (l + 1))
end

function _big_spectral_matrix(
    c::Complex{T},
    s::Int,
    m::Int,
    size::Int,
) where {T<:AbstractFloat}
    lmin = max(abs(m), abs(s))
    matrix = zeros(Complex{T}, size, size)
    @inbounds for l in lmin:(lmin + size - 1)
        row = l - lmin + 1
        for lprime in l:min(l + 2, lmin + size - 1)
            column = lprime - lmin + 1
            value = if lprime == l
                f = _big_swsh_f(T, s, l, m)
                g = _big_swsh_g(T, s, l, m)
                h = _big_swsh_h(T, s, l, m)
                b = f * _big_swsh_g(T, s, l + 1, m) +
                    g * _big_swsh_f(T, s, l - 1, m) + h^2
                T(l * (l + 1) - s * (s + 1)) -
                    c^2 * b + T(2s) * c * h
            elseif lprime == l + 1
                g = _big_swsh_g(T, s, lprime, m)
                e = g * (_big_swsh_h(T, s, lprime - 1, m) +
                    _big_swsh_h(T, s, lprime, m))
                -c^2 * e + T(2s) * c * g
            else
                -c^2 * _big_swsh_g(T, s, lprime, m) *
                    _big_swsh_g(T, s, lprime - 1, m)
            end
            matrix[row, column] = value
            matrix[column, row] = value
        end
    end
    return matrix
end

function _refine_angular_pair(
    seed::AngularEigenpair{Float64},
    precision_bits::Int;
    max_iterations::Int=20,
)
    precision_bits >= 64 ||
        throw(ArgumentError("precision_bits must be at least 64."))
    return setprecision(precision_bits) do
        T = BigFloat
        c = Complex{T}(T(real(seed.c)), T(imag(seed.c)))
        matrix = _big_spectral_matrix(
            c, seed.s, seed.m, seed.matrix_size)
        coefficients = Complex{T}.(seed.coefficients)
        pivot = argmax(abs.(coefficients))
        coefficients ./= coefficients[pivot]
        value = Complex{T}(
            T(real(seed.angular_sep)), T(imag(seed.angular_sep)))
        target = T(2)^(-T(precision_bits - 16))

        for _ in 1:max_iterations
            residual_vector = matrix * coefficients - value * coefficients
            normalization_error = coefficients[pivot] - one(Complex{T})
            residual = max(norm(residual_vector), abs(normalization_error))
            residual <= target && break

            size = seed.matrix_size
            jacobian = zeros(Complex{T}, size + 1, size + 1)
            jacobian[1:size, 1:size] .= matrix
            @inbounds for index in 1:size
                jacobian[index, index] -= value
                jacobian[index, size + 1] = -coefficients[index]
            end
            jacobian[size + 1, pivot] = one(Complex{T})
            rhs = -vcat(residual_vector, normalization_error)
            correction = jacobian \ rhs
            coefficients .+= view(correction, 1:size)
            value += correction[size + 1]
        end

        coefficients ./= norm(coefficients)
        reference = Complex{T}.(seed.coefficients)
        phase, overlap = _phase_align!(coefficients, reference)
        lambda = value + c^2 - T(2 * seed.m) * c
        residual = _relative_eigen_residual(
            matrix, value, coefficients)
        condition = _complex_symmetric_eigenvector_condition(coefficients)
        return AngularEigenpair(
            seed.s, seed.l, seed.m, c, seed.sheet_id, value, lambda,
            coefficients, seed.matrix_size, seed.truncation_order,
            precision_bits, :parallel_transport, residual, overlap,
            T(seed.overlap_margin), T(seed.spectral_gap), condition, phase,
            :high_precision_refined)
    end
end

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
    high_pair = continue_angular_mode(
        s, l, m, c;
        sheet_id=Symbol(sheet_id, :_order_high),
        truncation_order=high_order, cache=nothing, kwargs...)
    refined_low = _refine_angular_pair(high_pair, low_bits)
    refined_high = _refine_angular_pair(high_pair, high_bits)
    low_coefficients = _resize_coefficients(
        low_pair.coefficients, length(high_pair.coefficients))
    vector_overlap = abs(dot(low_coefficients, high_pair.coefficients)) /
        max(norm(low_coefficients) * norm(high_pair.coefficients),
            floatmin(Float64))
    truncation_drift = abs(low_pair.lambda - high_pair.lambda) /
        max(1.0, abs(low_pair.lambda), abs(high_pair.lambda))
    precision_drift = Float64(abs(refined_low.lambda - refined_high.lambda) /
        max(one(BigFloat), abs(refined_low.lambda), abs(refined_high.lambda)))
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

function spin_weighted_spheroidal_harmonic(
    pair::AngularEigenpair;
    method="auto",
)
    formatted_method = _format_method_name(method)
    formatted_method == "chebyshev" && throw(ArgumentError(
        "continued complex angular eigenpairs require a spectral harmonic backend."))
    coefficients_params = SpectralDecompositionInputParams(
        pair.s, pair.l, pair.m, pair.c, pair.matrix_size)
    l_list = construct_all_l_in_matrix(
        pair.s, pair.m, pair.matrix_size)
    spherical_harmonics_l = [
        abs(pair.coefficients[index]) >= _TOLERANCE ?
        spin_weighted_spherical_harmonic(
            pair.s, l_list[index], pair.m; method=formatted_method) :
        nothing for index in eachindex(l_list)
    ]
    return SpinWeightedSpheroidalHarmonicFunction(
        coefficients_params, pair.coefficients, spherical_harmonics_l,
        one(real(eltype(pair.coefficients))), pair.lambda, :spectral, nothing)
end

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
