# Spectral decomposition of the angular equation in spin-weighted spherical harmonics: the band
# matrix, its eigenpairs (found in Float64 with LAPACK and refined in a wider working precision), and
# the matrix size the coefficient tail asks for.
using LinearAlgebra
using LinearAlgebra.BLAS: @blasfunc, BlasInt

const SWSH_LAPACK = LinearAlgebra.LAPACK.liblapack
const SWSH_MIN_BUFFER = 10
const SWSH_MAX_REFINEMENTS = 20
# Number of trailing coefficients checked to decide whether the truncation is adequate
const SWSH_COEFFICIENT_TAIL_WINDOW = 4

# Below this |c| the eigenvalue is taken from its expansion to second order in c (_small_c_lambda)
# instead of a banded solve: in Float64 the neglected O(c^3) terms are then about 1e-18, below the
# resolution of an eigenvalue of order one
_small_c_limit(::Type{T}) where {T} = _tol(T, 1e-6)
# The eigenvalue-only solve accepts a matrix size once its eigenvalue agrees with the one at the
# previous size to within atol + rtol * max(|previous|, |current|): a few hundred eps absolute and a
# few tens of eps relative (in Float64), so that roundoff alone cannot keep the refinement going
_eigenvalue_atol(::Type{T}) where {T} = _tol(T, 5e-14)
_eigenvalue_rtol(::Type{T}) where {T} = _tol(T, 5e-15)
# Allowed deviation 1 - |overlap| of a computed eigenvector from the exact one; it enters only through
# _eigenvector_gap_tol below
_eigenvector_overlap_tol(::Type{T}) where {T} = _tol(T, 5e-13)
# Largest relative residual ||M v - lambda v|| / ((||M|| + |lambda|) ||v||) accepted for a banded eigenpair
# (see _real_lambda_residual); a converged LAPACK eigenpair sits at the eps level, far below
_eigenvector_residual_tol(::Type{T}) where {T} = _tol(T, 5e-13)
# Target size of the neglected spectral coefficients, used to estimate the matrix size
_coefficient_tail_tol(::Type{T}) where {T} = _tol(T, 1e-14)

#=
The truncation is adequate when the last few coefficients are below the coefficient tail
tolerance: the coefficients fall off faster than exponentially, so the neglected ones beyond N
are smaller still.
=#
function _coefficient_tail(coefficients)
    N = length(coefficients)
    return norm(view(coefficients, max(1, N - SWSH_COEFFICIENT_TAIL_WINDOW + 1):N))
end
#=
Roundoff perturbs an eigenvector by an angle of about eps*||M||/gap, where gap is the
distance to the nearest other eigenvalue. Below this relative gap that angle exceeds the
one allowed by the eigenvector overlap tolerance (1 - overlap ≈ angle^2/2), so the eigenvector
cannot be resolved in T even though the eigenvalue still can.
=#
_eigenvector_gap_tol(::Type{T}) where {T} = eps(T) / sqrt(2 * _eigenvector_overlap_tol(T))
# Working precision (bits) in which a Float64 eigenvalue is recomputed when its mode cannot be
# followed in Float64 (see Teukolsky_lambda_const)
const SWSH_NEAR_DEGENERATE_PRECISION = 128

#=
The eigenvalues and eigenvectors are first found in Float64 with LAPACK (the seed). In Float64
the seed is the result; a wider working precision refines it. The precision of the seed is that of
Float64, so a type no wider than that has nothing to refine.
=#
_seed_is_final(::Type{T}) where {T} = precision(T) <= precision(Float64)

_precision_name(::Type{Float64}) = "Float64"
_precision_name(::Type{T}) where {T} = "$(precision(T))-bit $T"

#=
The eigenvalue is accurate to about eps*||M||, so while the gap (relative to ||M||) is
comfortably larger than that, the eigenvalue still tells the two modes apart.
It can then seed Leaver's method, which resolves the harmonic in a wider precision.
Unseeded, Leaver's march in c may give up here, in any precision, once the neighbouring mode comes
too close to follow (see _leaver_march); the seed lets the solve skip the march.
=#
function _unresolvable_eigenvector_message(s, l, m, c, relative_gap)
    T = real(typeof(float(c)))
    message = "The SWSH eigenvector for (s, l, m) = ($s, $l, $m), c = $c cannot be " *
        "resolved in $(_precision_name(T)): the mode is nearly degenerate with a neighbouring mode " *
        "(relative gap $relative_gap, below $(_eigenvector_gap_tol(T))), so the " *
        "harmonic would be an arbitrary mix of the two. "
    if relative_gap > 10eps(T)
        return message * "Pass c in a wider precision (for example precision=256), which resolves the " *
            "mode. The eigenvalue is still accurate, see " *
            "spin_weighted_spheroidal_eigenvalue; the harmonic can also be refined in a wider " *
            "precision with Leaver's method, seeded with that eigenvalue: " *
            "spin_weighted_spheroidal_harmonic($s, $l, $m, big($c); method=\"leaver\", " *
            "lambda0=big(spin_weighted_spheroidal_eigenvalue($s, $l, $m, $c)))."
    end
    return message * "Here even the eigenvalue cannot tell the two modes apart in " *
        "$(_precision_name(T)), so the mode cannot be identified without a wider precision throughout: " *
        "pass c in a wider precision, for example precision=256."
end

function _determine_matrix_size_N(s::Int, l::Int, m::Int, c=0)
    #=
    Estimate a suitable value of N for the spectral decomposition

    N covers every spherical harmonic up to the target mode, plus the harmonics
    beyond it that the mode still needs, plus a 'buffer' of SWSH_MIN_BUFFER (in Float64).
    At large |c| the mode is confined near a pole with a width ~ 1/sqrt(|c|), so its
    spherical-harmonic coefficients fall off like a Gaussian in l of width ~ sqrt(|c|):
    about sqrt(2 log(1/tol) |c|) more harmonics bring them below tol. This is an
    estimate; the real-c solvers that start from it check the coefficient tail afterwards and
    enlarge N if needed, while the complex-c continuation uses it as is.
    =#
    T = _float_type(c)
    N = l - max(abs(m), abs(s)) + 1
    beyond = ceil(Int, sqrt(2log(1 / _coefficient_tail_tol(T)) * abs(c)))
    # the buffer covers a geometric fall-off, so it grows with the number of digits
    return N + beyond + _nterms(T, SWSH_MIN_BUFFER)
end

function Fslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    # 'Edge' case where l is -1, this can happen when both |m| and |s| are 0 (since lmin = max(|m|, |s|))
    (l == -1 && abs(m) == 0 && abs(s) == 0) && return zero(T)
    lp1 = T(l + 1)
    return sqrt(((lp1^2 - T(m)^2) / (T(2l + 3) * T(2l + 1))) *
        ((lp1^2 - T(s)^2) / lp1^2))
end

function Gslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    iszero(l) && return zero(T)
    ll = T(l)
    return sqrt(((ll^2 - T(m)^2) / (4ll^2 - 1)) *
        ((ll^2 - T(s)^2) / ll^2))
end

function Hslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    (l == 0 || s == 0) ? zero(T) : T(-m*s)/T(l*(l+1))
end

function Aslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    Fslm(T, s, l, m)*Fslm(T, s, l+1, m)
end

function Bslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    Fslm(T, s, l, m)*Gslm(T, s, l+1, m) + Gslm(T, s, l, m)*Fslm(T, s, l-1, m) + Hslm(T, s, l, m)^2
end

function Cslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    Gslm(T, s, l, m)*Gslm(T, s, l-1, m)
end

function Dslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    Fslm(T, s, l, m)*(Hslm(T, s, l+1, m) + Hslm(T, s, l, m))
end

function Eslm(::Type{T}, s::Int, l::Int, m::Int) where {T}
    Gslm(T, s, l, m)*(Hslm(T, s, l-1, m) + Hslm(T, s, l, m))
end

function eigenvalue_Schwarzschild(s::Int, l::Int)
    l*(l+1) - s*(s+1)
end

@inline function _small_c_lambda(c::T, s::Int, l::Int, m::Int) where {T<:Real}
    g_lower = Gslm(T, s, l, m)
    g_upper = Gslm(T, s, l + 1, m)
    h = Hslm(T, s, l, m)
    lower_shift = iszero(l) ? zero(T) : g_lower^2 / l
    lambda0 = T(eigenvalue_Schwarzschild(s, l))
    lambda1 = 2s * h - 2m
    lambda2 = 1 - Bslm(T, s, l, m) +
        2s^2 * (lower_shift - g_upper^2 / (l + 1))
    return muladd(c, muladd(c, lambda2, lambda1), lambda0)
end

function spectral_matrix_coefficient(c, s::Int, m::Int, l::Int, lprime::Int)
    # This is Eq (55)
    T = _float_type(c)
    if lprime == l-2
        return -c^2 * Aslm(T, s, lprime, m)
    elseif lprime == l-1
        return -c^2 * Dslm(T, s, lprime, m) + 2*c*s*Fslm(T, s, lprime, m)
    elseif lprime == l
        return eigenvalue_Schwarzschild(s, lprime) - c^2 * Bslm(T, s, lprime, m) + 2*c*s*Hslm(T, s, lprime, m)
    elseif lprime == l+1
        return -c^2 * Eslm(T, s, lprime, m) + 2*c*s*Gslm(T, s, lprime, m)
    elseif lprime == l+2
        return -c^2 * Cslm(T, s, lprime, m)
    else
        return 0
    end
end

function construct_all_l_in_matrix(s::Int, m::Int, N::Int)
    lmin = max(abs(m), abs(s))
    lmax = N + lmin - 1
    [k for k in lmin:lmax]
end

function construct_spectral_matrix(c, s::Int, m::Int, N::Int)
    all_l_in_matrix = construct_all_l_in_matrix(s, m, N)
    lmin = all_l_in_matrix[1]
    lmax = all_l_in_matrix[N]

    # Matrix is symmetric, fill in only the upper right portion
    spectral_matrix = zeros(Complex{_float_type(c)}, N, N)
    for i in lmin:lmax
        for j in i:lmax
            spectral_matrix[i-lmin+1, j-lmin+1] = spectral_matrix_coefficient(c, s, m, i, j)
        end
    end
    spectral_matrix = Array(Symmetric(spectral_matrix))
end

# The spectral matrix in upper band storage, band[3 + i - j, j] = M[i, j], for real or complex c
function _spectral_band(c, s::Int, m::Int, N::Int)
    T = _float_type(c)
    lmin = max(abs(m), abs(s))
    lmax = lmin + N - 1
    bandwidth = 2
    band = zeros(typeof(c * one(T)), bandwidth + 1, N)
    @inbounds for l in lmin:lmax
        i = l - lmin + 1
        for lprime in l:min(l + bandwidth, lmax)
            j = lprime - lmin + 1
            value = if lprime == l
                spherical = T(eigenvalue_Schwarzschild(s, l))
                muladd(c^2, 1 - Bslm(T, s, l, m),
                    muladd(2c, s * Hslm(T, s, l, m) - m, spherical))
            else
                spectral_matrix_coefficient(c, s, m, l, lprime)
            end
            band[bandwidth + 1 + i - j, j] = value
        end
    end
    return band
end

# The real band matrix of real c (the λ-matrix, M + c^2 - 2mc), which LAPACK's dsbevx takes
_real_lambda_band(c::Real, s::Int, m::Int, N::Int) = _spectral_band(float(c), s, m, N)

_selected_banded_eigenvalue!(band::Matrix{Float64}, index::Int) =
    first(_selected_banded_eigenvalues!(band, index, index))

# Eigenvalues lower:upper (in ascending order) of the symmetric band matrix, without eigenvectors
function _selected_banded_eigenvalues!(band::Matrix{Float64}, lower::Int, upper::Int)
    n = BlasInt(size(band, 2))
    kd = BlasInt(size(band, 1) - 1)
    leading_band = BlasInt(size(band, 1))
    q = zeros(Float64, 1)
    leading_q = BlasInt(1)
    lower_value = 0.0
    upper_value = 0.0
    lower_index = BlasInt(lower)
    upper_index = BlasInt(upper)
    expected = upper - lower + 1
    # LAPACK recommends twice the safe minimum for high relative accuracy.
    absolute_tolerance = 2floatmin(Float64)
    found = Ref{BlasInt}()
    eigenvalues = zeros(Float64, Int(n))
    eigenvectors = zeros(Float64, 1)
    leading_eigenvectors = BlasInt(1)
    work = zeros(Float64, 7 * Int(n))
    integer_work = zeros(BlasInt, 5 * Int(n))
    failed = zeros(BlasInt, Int(n))
    info = Ref{BlasInt}()

    ccall((@blasfunc(dsbevx_), SWSH_LAPACK), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Clong, Clong, Clong),
        'N', 'I', 'U', n, kd, band, leading_band, q, leading_q,
        lower_value, upper_value, lower_index, upper_index,
        absolute_tolerance, found, eigenvalues, eigenvectors,
        leading_eigenvectors, work, integer_work, failed, info, 1, 1, 1)

    info[] == 0 || error("LAPACK dsbevx failed with info=$(info[]).")
    found[] == expected ||
        error("LAPACK dsbevx returned $(found[]) eigenvalues instead of $expected.")
    return eigenvalues[1:expected]
end

# Infinity norm of the symmetric matrix held in upper band storage
function _symmetric_band_inf_norm(band::AbstractMatrix)
    kd = size(band, 1) - 1
    n = size(band, 2)
    row_sums = zeros(real(eltype(band)), n)
    @inbounds for j in 1:n, i in max(1, j - kd):j
        value = abs(band[kd + 1 + i - j, j])
        row_sums[i] += value
        i != j && (row_sums[j] += value)
    end
    return maximum(row_sums)
end

#=
Return the requested eigenpair and, with neighbours, the gap to the nearest neighbouring
eigenvalue relative to the matrix norm (NaN otherwise). The neighbours come from the same
LAPACK call, which costs little more than the eigenpair alone: most of the cost of a call
is the reduction to tridiagonal form, which they share.
=#
function _selected_banded_eigenpair!(band::Matrix{Float64}, index::Int;
        neighbours::Bool=false)
    n = BlasInt(size(band, 2))
    kd = BlasInt(size(band, 1) - 1)
    leading_band = BlasInt(size(band, 1))
    matrix_norm = neighbours ? _symmetric_band_inf_norm(band) : NaN # dsbevx overwrites band
    q = zeros(Float64, Int(n), Int(n))
    leading_q = n
    lower_value = 0.0
    upper_value = 0.0
    lower_index = BlasInt(neighbours ? max(1, index - 1) : index)
    upper_index = BlasInt(neighbours ? min(Int(n), index + 1) : index)
    expected = Int(upper_index - lower_index) + 1
    absolute_tolerance = 2floatmin(Float64)
    found = Ref{BlasInt}()
    eigenvalues = zeros(Float64, Int(n))
    eigenvectors = zeros(Float64, Int(n), expected)
    leading_eigenvectors = n
    work = zeros(Float64, 7 * Int(n))
    integer_work = zeros(BlasInt, 5 * Int(n))
    failed = zeros(BlasInt, Int(n))
    info = Ref{BlasInt}()

    ccall((@blasfunc(dsbevx_), SWSH_LAPACK), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{Float64}, Ref{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ref{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt},
            Ref{BlasInt}, Clong, Clong, Clong),
        'V', 'I', 'U', n, kd, band, leading_band, q, leading_q,
        lower_value, upper_value, lower_index, upper_index,
        absolute_tolerance, found, eigenvalues, eigenvectors,
        leading_eigenvectors, work, integer_work, failed, info, 1, 1, 1)

    info[] == 0 || error("LAPACK dsbevx failed with info=$(info[]).")
    found[] == expected ||
        error("LAPACK dsbevx returned $(found[]) eigenpairs instead of $expected.")
    position = index - Int(lower_index) + 1 # the requested mode within the block
    value = eigenvalues[position]
    relative_gap = if neighbours
        gap = minimum((abs(eigenvalues[k] - value) for k in 1:expected
            if k != position); init=Inf)
        gap / max(matrix_norm, floatmin(Float64))
    else
        NaN
    end
    vector = view(eigenvectors, :, position)
    vector ./= norm(vector)
    pivot = vector[index]
    if iszero(pivot)
        pivot = vector[argmax(abs.(vector))]
    end
    signbit(pivot) && (vector .*= -1)
    return value, ComplexF64.(vector), relative_gap
end

#=
Linear algebra in the working precision on the band matrix: LU with partial pivoting of a
shifted band matrix (O(N) per solve, as LAPACK's gbtrf/gbtrs), Rayleigh quotient iteration for
an eigenpair, and Sturm counts (the inertia of the LDLᵀ factorization) to select an eigenvalue by
its position.
=#

# LU factors of M - shift*I for M in symmetric upper band storage (kd super-diagonals), in the
# general band layout of gbtrf: A[i, j] is held at lu[2kd + 1 + i - j, j]
function _shifted_band_lu(band::AbstractMatrix, shift)
    kd = size(band, 1) - 1
    n = size(band, 2)
    S = promote_type(eltype(band), typeof(shift))
    lu = zeros(S, 3kd + 1, n)
    @inbounds for j in 1:n, i in max(1, j - kd):min(n, j + kd)
        value = i <= j ? band[kd + 1 + i - j, j] : band[kd + 1 + j - i, i]
        i == j && (value -= shift)
        lu[2kd + 1 + i - j, j] = value
    end
    pivots = zeros(Int, n)
    kv = 2kd
    @inbounds for j in 1:n
        km = min(kd, n - j)
        p = 0
        largest = -one(real(S))
        for t in 0:km
            size_t = abs(lu[kv + 1 + t, j])
            size_t > largest && (largest = size_t; p = t)
        end
        pivots[j] = j + p
        last_column = min(j + kv, n)
        if p != 0
            for col in j:last_column
                a, b = kv + 1 + j - col, kv + 1 + j + p - col
                lu[a, col], lu[b, col] = lu[b, col], lu[a, col]
            end
        end
        pivot = lu[kv + 1, j]
        iszero(pivot) && (lu[kv + 1, j] = pivot = eps(real(S)) * max(largest, floatmin(real(S))))
        for t in 1:km
            lu[kv + 1 + t, j] /= pivot
        end
        for col in j+1:last_column, t in 1:km
            lu[kv + 1 + j + t - col, col] -= lu[kv + 1 + t, j] * lu[kv + 1 + j - col, col]
        end
    end
    return lu, pivots, kd
end

function _shifted_band_solve((lu, pivots, kd), rhs::AbstractVector)
    n = length(rhs)
    kv = 2kd
    x = Vector{promote_type(eltype(lu), eltype(rhs))}(rhs)
    @inbounds for j in 1:n
        p = pivots[j]
        p != j && ((x[j], x[p]) = (x[p], x[j]))
        for t in 1:min(kd, n - j)
            x[j + t] -= lu[kv + 1 + t, j] * x[j]
        end
    end
    @inbounds for j in n:-1:1
        x[j] /= lu[kv + 1, j]
        for i in max(1, j - kv):j-1
            x[i] -= lu[kv + 1 + i - j, j] * x[j]
        end
    end
    return x
end

# M * v for M in symmetric upper band storage (unconjugated, so also for complex symmetric M)
function _symmetric_band_mul(band::AbstractMatrix, v::AbstractVector)
    kd = size(band, 1) - 1
    n = size(band, 2)
    product = zeros(promote_type(eltype(band), eltype(v)), n)
    @inbounds for j in 1:n, i in max(1, j - kd):j
        value = band[kd + 1 + i - j, j]
        product[i] += value * v[j]
        i != j && (product[j] += value * v[i])
    end
    return product
end

function _band_relative_residual(band, value, vector)
    residual = norm(_symmetric_band_mul(band, vector) - value * vector)
    scale = (_symmetric_band_inf_norm(band) + abs(value)) * norm(vector)
    return residual / max(scale, floatmin(real(eltype(band))))
end

#=
Rayleigh quotient iteration from an accurate seed: the shifted solve is nearly singular by
design, and each step multiplies the number of correct digits by about three (two for a complex
symmetric matrix, whose quotient vᵀMv/vᵀv is not conjugated). It stops once the residual has
reached roundoff in the working precision or stops improving.
=#
function _rayleigh_quotient_iteration(band::AbstractMatrix, value, vector; max_iterations::Int=12)
    R = real(eltype(band))
    target = 16 * sqrt(R(length(vector))) * eps(R)
    vector = vector / norm(vector)
    residual = _band_relative_residual(band, value, vector)
    for _ in 1:max_iterations
        residual <= target && break
        next = _shifted_band_solve(_shifted_band_lu(band, value), vector)
        next ./= norm(next)
        next_value = transpose(next) * _symmetric_band_mul(band, next) / (transpose(next) * next)
        next_value = eltype(band) <: Real ? real(next_value) : next_value
        next_residual = _band_relative_residual(band, next_value, next)
        next_residual < residual || break
        vector, value, residual = next, next_value, next_residual
    end
    return value, vector, residual
end

# Number of eigenvalues of the real symmetric band matrix below x: the negative pivots of the
# LDLᵀ factorization of M - x I (Sylvester's law of inertia)
function _sturm_count(band::AbstractMatrix{R}, x::R) where {R<:Real}
    kd = size(band, 1) - 1
    n = size(band, 2)
    tiny = eps(R) * max(_symmetric_band_inf_norm(band), floatmin(R))
    # rows of L (below the diagonal) and the pivots, built column by column
    L = zeros(R, kd, n)   # L[t, j] = l[j + t, j]
    d = zeros(R, n)
    count = 0
    @inbounds for j in 1:n
        dj = band[kd + 1, j] - x
        for k in max(1, j - kd):j-1
            dj -= L[j - k, k]^2 * d[k]
        end
        iszero(dj) && (dj = -tiny)
        d[j] = dj
        dj < 0 && (count += 1)
        for t in 1:min(kd, n - j)
            i = j + t
            value = band[kd + 1 + j - i, i]
            for k in max(1, i - kd):j-1
                value -= L[i - k, k] * L[j - k, k] * d[k]
            end
            L[t, j] = value / dj
        end
    end
    return count
end

# The index-th eigenvalue (ascending) of the real symmetric band matrix, by bisection from a
# bracket around a seed, to roundoff in the working precision
function _sturm_eigenvalue(band::AbstractMatrix{R}, index::Int, seed) where {R<:Real}
    scale = max(_symmetric_band_inf_norm(band), floatmin(R))
    width = 4eps(Float64) * scale
    lower, upper = R(seed) - width, R(seed) + width
    while _sturm_count(band, lower) >= index
        lower -= width; width *= 2
    end
    while _sturm_count(band, upper) < index
        upper += width; width *= 2
    end
    while upper - lower > 2eps(R) * max(abs(lower), abs(upper), scale)
        middle = (lower + upper) / 2
        (middle == lower || middle == upper) && break
        _sturm_count(band, middle) >= index ? (upper = middle) : (lower = middle)
    end
    return (lower + upper) / 2
end

# Unit vector with the target coefficient positive (or the largest one, if that vanishes)
function _sign_normalized(vector::AbstractVector, index::Int)
    vector = vector / norm(vector)
    pivot = vector[index]
    iszero(pivot) && (pivot = vector[argmax(abs.(vector))])
    signbit(real(pivot)) && (vector = -vector)
    return vector
end

_padded(vector::AbstractVector{S}, N::Int) where {S} =
    N <= length(vector) ? vector[1:N] : vcat(vector, zeros(S, N - length(vector)))

#=
The eigenpair at index `index` of the real band matrix in the working precision, from the Float64
seed (λ, v). When the seed identifies the mode, Rayleigh quotient iteration refines it. When it does
not (the gap to a neighbour is below what Float64 resolves), the eigenvalue is selected by its
position with Sturm counts, and the eigenvector follows from inverse iteration; the neighbours are
located the same way to measure the gap, which must then be resolvable in the working precision.
Return the eigenvalue, the unit eigenvector with the target coefficient positive, and the relative gap.
=#
function _refined_real_eigenpair(band::AbstractMatrix{R}, index::Int, seed_value, seed_vector,
        seed_gap, identified::Bool, s, l, m, c) where {R<:Real}
    if identified
        value, vector, _ = _rayleigh_quotient_iteration(band, R(seed_value), R.(real.(seed_vector)))
        return value, _sign_normalized(vector, index), R(seed_gap)
    end
    n = size(band, 2)
    value = _sturm_eigenvalue(band, index, seed_value)
    scale = max(_symmetric_band_inf_norm(band), floatmin(R))
    neighbours = [_sturm_eigenvalue(band, k, seed_value) for k in (index - 1, index + 1) if 1 <= k <= n]
    gap = minimum((abs(v - value) for v in neighbours); init=R(Inf)) / scale
    gap >= _eigenvector_gap_tol(R) || error(_unresolvable_eigenvector_message(s, l, m, c, gap))
    vector = R.(real.(seed_vector))
    factors = _shifted_band_lu(band, value)
    for _ in 1:3
        vector = _shifted_band_solve(factors, vector)
        vector ./= norm(vector)
    end
    value, vector, _ = _rayleigh_quotient_iteration(band, value, vector)
    return value, _sign_normalized(vector, index), gap
end

#=
Also return the relative gap to the neighbouring eigenvalues. With check_gap, throw if
the eigenvector cannot be resolved (see _eigenvector_gap_tol).
=#
function _real_lambda_eigenpair_at_size(
    c::Real,
    s::Int,
    l::Int,
    m::Int,
    N::Int;
    check_gap::Bool=true,
)
    T = _float_type(c)
    index = _ell_index_in_matrix(s, l, m, N)
    lambda, coefficients, relative_gap = _selected_banded_eigenpair!(
        _real_lambda_band(Float64(c), s, m, N), index; neighbours=true)
    if _seed_is_final(T)
        check_gap && relative_gap < _eigenvector_gap_tol(Float64) &&
            error(_unresolvable_eigenvector_message(s, l, m, c, relative_gap))
        return lambda, coefficients, relative_gap
    end
    identified = relative_gap >= _eigenvector_gap_tol(Float64)
    value, vector, gap = _refined_real_eigenpair(_real_lambda_band(T(c), s, m, N), index,
        lambda, coefficients, relative_gap, identified, s, l, m, c)
    return value, Complex{T}.(vector), gap
end

function _real_lambda_residual(c::Real, s::Int, m::Int, lambda,
        coefficients)
    T = _float_type(c)
    N = length(coefficients)
    lmin = max(abs(m), abs(s))
    lmax = lmin + N - 1
    product = zeros(Complex{T}, N)
    maximum_row_sum = zero(T)
    c = T(c)
    @inbounds for l in lmin:lmax
        i = l - lmin + 1
        row_sum = zero(T)
        for lprime in l:min(l + 2, lmax)
            j = lprime - lmin + 1
            value = if lprime == l
                spherical = T(eigenvalue_Schwarzschild(s, l))
                muladd(c^2, 1 - Bslm(T, s, l, m),
                    muladd(2c, s * Hslm(T, s, l, m) - m, spherical))
            else
                real(spectral_matrix_coefficient(c, s, m, l, lprime))
            end
            product[i] += value * coefficients[j]
            row_sum += abs(value)
            if j != i
                product[j] += value * coefficients[i]
            end
        end
        maximum_row_sum = max(maximum_row_sum, row_sum)
    end
    residual = norm(product - lambda * coefficients)
    scale = (maximum_row_sum + abs(lambda)) * norm(coefficients)
    return residual / max(scale, floatmin(T))
end

#=
Start from the size estimated from |c|, which is usually enough, and accept the first
solve whose coefficient tail and residual are small: once the neglected coefficients are
below tolerance, a larger matrix could only change the eigenpair at that level. Increase
the size otherwise, whether it is the tail or the residual that is still too large. The seed does not throw: whether its
eigenvector is resolved is decided by the caller, in the working precision.
=#
function _real_eigenpair_seed(c::Float64, s::Int, l::Int, m::Int)
    if iszero(c)
        size = _determine_matrix_size_N(s, l, m)
        index = _ell_index_in_matrix(s, l, m, size)
        coefficients = zeros(ComplexF64, size)
        coefficients[index] = 1.0 + 0.0im
        return (lambda=Float64(eigenvalue_Schwarzschild(s, l)),
            coefficients, size, refinement=0, tail=0.0, residual=0.0, gap=Inf, converged=true)
    end

    size = _determine_matrix_size_N(s, l, m, c)
    # Each refinement adds at least 8 harmonics, and |c|/8 once |c| exceeds 64, so that few solves
    # are needed even where the estimate falls well short
    step = max(8, ceil(Int, abs(c) / 8))
    for refinement in 0:SWSH_MAX_REFINEMENTS
        lambda, coefficients, gap = _real_lambda_eigenpair_at_size(
            c, s, l, m, size; check_gap=false)
        tail = _coefficient_tail(coefficients)
        residual = _real_lambda_residual(c, s, m, lambda, coefficients)
        converged = tail <= _coefficient_tail_tol(Float64) &&
            residual <= _eigenvector_residual_tol(Float64)
        (converged || refinement == SWSH_MAX_REFINEMENTS) &&
            return (; lambda, coefficients, size, refinement, tail, residual, gap, converged)
        size += step
    end
end

function _adaptive_real_eigenpair(c::Real, s::Int, l::Int, m::Int)
    T = _float_type(c)
    seed = _real_eigenpair_seed(Float64(c), s, l, m)
    identified = seed.gap >= _eigenvector_gap_tol(Float64)
    if _seed_is_final(T)
        # An unresolvable eigenvector also shows up as one that never settles
        identified || error(_unresolvable_eigenvector_message(s, l, m, c, seed.gap))
        seed.converged || error("SWSH eigenpair failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
        return (; seed.lambda, seed.coefficients, seed.size, seed.refinement, seed.tail, seed.residual, seed.gap)
    end
    identified && !seed.converged &&
        error("SWSH eigenpair failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
    c = T(c)
    index = _ell_index_in_matrix(s, l, m, seed.size)
    size = max(seed.size, _determine_matrix_size_N(s, l, m, c))
    step = max(8, ceil(Int, abs(c) / 8)) # as in _real_eigenpair_seed
    lambda, vector = seed.lambda, seed.coefficients
    for refinement in 0:SWSH_MAX_REFINEMENTS
        lambda, vector, gap = _refined_real_eigenpair(_real_lambda_band(c, s, m, size), index,
            lambda, _padded(vector, size), seed.gap, identified, s, l, m, c)
        coefficients = Complex{T}.(vector)
        tail = _coefficient_tail(coefficients)
        if tail <= _coefficient_tail_tol(T)
            residual = _real_lambda_residual(c, s, m, lambda, coefficients)
            return (; lambda, coefficients, size, refinement, tail, residual, gap)
        end
        size += step
    end
    error("SWSH eigenpair failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

function _real_lambda_eigenvalue_at_size(
    c::Real,
    s::Int,
    l::Int,
    m::Int,
    N::Int,
)
    T = _float_type(c)
    index = _ell_index_in_matrix(s, l, m, N)
    seed = _selected_banded_eigenvalue!(_real_lambda_band(Float64(c), s, m, N), index)
    _seed_is_final(T) && return seed
    return _sturm_eigenvalue(_real_lambda_band(T(c), s, m, N), index, seed)
end

function _ell_index_in_matrix(s::Int, l::Int, m::Int, N::Int)
    lmin = max(abs(m), abs(s))
    idx = l - lmin + 1
    if idx < 1 || idx > N
        error("Target l=$l is outside the spectral matrix range for N=$N")
    end
    return idx
end

# Harmonics kept beyond the target mode when continuing to c; N = -1 estimates it from c
@inline function _complex_truncation_order(
    s::Int,
    l::Int,
    m::Int,
    N::Int,
    c,
)
    N == -1 && (N = _determine_matrix_size_N(s, l, m, c))
    order = N - (l - max(abs(m), abs(s)) + 1)
    order >= 0 ||
        throw(ArgumentError("N does not contain the target angular mode."))
    return order
end

function _spectral_backend(backend)
    value = Symbol(lowercase(String(backend)))
    value in (:auto, :banded, :dense) ||
        throw(ArgumentError(
            "backend must be :auto, :banded, or :dense."))
    return value
end

# Like _spectral_backend, but also rejects a backend that cannot handle c,
# instead of silently running a different solver
function _resolve_spectral_backend(backend, c)
    value = _spectral_backend(backend)
    if value == :banded && c isa Complex && !iszero(imag(c))
        throw(ArgumentError(
            "backend=:banded requires real c; use :auto or :dense."))
    end
    return value
end

#=
N = -1 picks the default size here, where it is known which branch is taken:
for complex c the same truncation as the default (:auto) route, so that the two
backends can be compared directly; for real c an estimate from |c| that is then
checked against the coefficient tail, and increased if needed.
=#
function _dense_spectral_decomposition(c, s::Int, l::Int, m::Int, N::Int=-1)
    if c isa Complex && !iszero(imag(c))
        # Complex eigenvalue ordering does not preserve the spherical mode label.
        pair = continue_angular_mode(
            s, l, m, float(c); backend=:dense,
            truncation_order=_complex_truncation_order(s, l, m, N, c))
        return pair.angular_sep, _overlap_phase_fixed(pair.coefficients, s, l, m)
    end
    N == -1 || return _dense_real_eigenpair_at_size(c, s, l, m, N)
    N = _determine_matrix_size_N(s, l, m, c)
    step = max(8, ceil(Int, abs(c) / 8)) # as in _real_eigenpair_seed
    for _ in 0:SWSH_MAX_REFINEMENTS
        angular_sep, coefficients = _dense_real_eigenpair_at_size(c, s, l, m, N)
        _coefficient_tail(coefficients) <= _coefficient_tail_tol(_float_type(c)) &&
            return angular_sep, coefficients
        N += step
    end
    error("SWSH dense eigenpair did not reach a small enough coefficient tail after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

#=
The dense backend solves the full eigenproblem in Float64 (for real c the mode is the
index-th eigenvalue); a wider working precision refines that eigenpair on the band matrix
exactly as the banded backend does, so the two backends differ only in the seed.
=#
function _dense_real_eigenpair_at_size(c, s::Int, l::Int, m::Int, N::Int)
    T = _float_type(c)
    idx = _ell_index_in_matrix(s, l, m, N)
    if c == 0
        coefficients = zeros(Complex{T}, N)
        coefficients[idx] = one(T) + zero(T)im
        return T(eigenvalue_Schwarzschild(s, l)), coefficients
    end
    seed_c = c isa Complex ? ComplexF64(c) : Float64(c)
    spectral_matrix = construct_spectral_matrix(seed_c, s, m, N)
    decomposition = eigen(spectral_matrix)
    relative_gap = _relative_spectral_gap(
        decomposition.values, idx, opnorm(spectral_matrix, Inf))
    angular_sep = decomposition.values[idx]
    vector = decomposition.vectors[:, idx]
    if _seed_is_final(T)
        relative_gap >= _eigenvector_gap_tol(Float64) ||
            error(_unresolvable_eigenvector_message(s, l, m, c, relative_gap))
        pivot = vector[idx]
        !iszero(pivot) && (vector /= pivot)
        coefficients = vector / sqrt(dot(vector, vector))
        return angular_sep, coefficients
    end
    # The band matrix holds the λ-matrix (angular_sep + c^2 - 2mc); the dense one the angular_sep matrix
    c, c64 = T(real(c)), Float64(real(c))
    shift = muladd(c, c, -2m * c)
    lambda, refined, _ = _refined_real_eigenpair(_real_lambda_band(c, s, m, N), idx,
        real(angular_sep) + muladd(c64, c64, -2m * c64), real.(vector), relative_gap,
        relative_gap >= _eigenvector_gap_tol(Float64), s, l, m, c)
    return Complex{T}(lambda - shift), Complex{T}.(refined)
end

function _spectral_decomposition(c, s::Int, l::Int, m::Int, N::Int=-1;
        backend=:auto)
    selected_backend = _resolve_spectral_backend(backend, c)
    if selected_backend != :dense && isreal(c)
        lambda, coefficients = if N == -1
            pair = _adaptive_real_eigenpair(real(c), s, l, m)
            pair.lambda, pair.coefficients
        else
            _real_lambda_eigenpair_at_size(real(c), s, l, m, N)
        end
        T = _float_type(c)
        angular_sep = lambda -
            muladd(T(real(c)), T(real(c)),
                -2m * T(real(c)))
        return angular_sep, coefficients
    end
    return _dense_spectral_decomposition(c, s, l, m, N)
end

function _adaptive_real_lambda(c::Real, s::Int, l::Int, m::Int)
    T = _float_type(c)
    iszero(c) && return (value=T(eigenvalue_Schwarzschild(s, l)),
        size=0, refinement=0, delta=zero(T))
    if abs(c) <= _small_c_limit(T)
        return (value=_small_c_lambda(T(c), s, l, m),
            size=_determine_matrix_size_N(s, l, m), refinement=0, delta=zero(T))
    end
    c = T(c)
    lmin = max(abs(m), abs(s))
    idx = l - lmin + 1
    # Starting size for the eigenvalue-only solve: ceil(|c|) + 8 harmonics beyond the target index, and
    # at least SWSH_MIN_BUFFER; the size is then enlarged (as in _real_eigenpair_seed) until the
    # eigenvalue settles
    size = max(idx + SWSH_MIN_BUFFER, idx + ceil(Int, abs(c)) + 8)
    step = max(8, ceil(Int, abs(c) / 8))
    previous = _real_lambda_eigenvalue_at_size(c, s, l, m, size)
    for refinement in 1:SWSH_MAX_REFINEMENTS
        next_size = size + step
        current = _real_lambda_eigenvalue_at_size(
            c, s, l, m, next_size)
        delta = abs(current - previous)
        threshold = _eigenvalue_atol(T) +
            _eigenvalue_rtol(T) * max(abs(previous), abs(current))
        delta <= threshold && return (value=current, size=next_size,
            refinement=refinement, delta=delta)
        size = next_size
        previous = current
    end
    error("SWSH eigenvalue failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

function _angular_eigenvalue(c::Real, s::Int, l::Int, m::Int, N::Int=-1)
    T = _float_type(c)
    iszero(c) && return T(eigenvalue_Schwarzschild(s, l))
    lambda = N == -1 ? _adaptive_real_lambda(c, s, l, m).value :
        _real_lambda_eigenvalue_at_size(c, s, l, m, N)
    return lambda - muladd(T(c), T(c), -2m * T(c))
end

function _angular_eigenvalue(c, s::Int, l::Int, m::Int, N::Int=-1)
    # Only the eigenvalue is needed, so take the eigenvalue-only real route,
    # which also works where the eigenvector cannot be resolved
    isreal(c) && return _angular_eigenvalue(real(c), s, l, m, N)
    if c isa Complex && !iszero(imag(c))
        return continue_angular_mode(
            s, l, m, c;
            truncation_order=_complex_truncation_order(s, l, m, N, c)).angular_sep
    end
    angular_sep, _ = _spectral_decomposition(c, s, l, m, N)
    return angular_sep
end

function angular_sep_const(c, s::Int, l::Int, m::Int, N::Int=-1)
    return _angular_eigenvalue(c, s, l, m, N)
end

function spectral_coefficients(c, s::Int, l::Int, m::Int, N::Int=-1;
        backend=:auto)
    selected_backend = _resolve_spectral_backend(backend, c)
    if c isa Complex && !iszero(imag(c)) &&
            selected_backend != :dense
        return _overlap_phase_fixed(continue_angular_mode(
            s, l, m, c;
            truncation_order=_complex_truncation_order(s, l, m, N, c)).coefficients, s, l, m)
    end
    _, coeffs = _spectral_decomposition(
        c, s, l, m, N; backend=selected_backend)
    return coeffs
end

@doc raw"""
    Teukolsky_lambda_const(c, s::Int, l::Int, m::Int, N::Int=-1)

Compute the eigenvalue ``\lambda`` of the spin-weighted spheroidal harmonic with spin weight `s`,
harmonic index `l`, azimuthal index `m` and spheroidicity `c` (``c = a\omega``) by spectral
decomposition, in the floating-point type of `c`. The convention is that of Teukolsky's radial
equation, ``\lambda = A + c^2 - 2mc`` with ``A`` the separation constant of the angular equation, so
that ``\lambda = l(l+1) - s(s+1)`` at `c = 0`.

`N` is the number of spin-weighted spherical harmonics in the decomposition; the default `N=-1`
determines a suitable value automatically. For real `c` only the eigenvalue is computed, so it is
available even where the eigenvector cannot be resolved. For complex `c`, "mode `l`" is the mode that
connects continuously to the spherical harmonic `l` along the straight path from `c = 0` (see
`continue_angular_mode`).

Note that `c` comes first, unlike in `spin_weighted_spheroidal_eigenvalue`, which is the
recommended interface; this function is kept for backward compatibility.
"""
function Teukolsky_lambda_const(c::Real, s::Int, l::Int, m::Int, N::Int=-1)
    T = _float_type(c)
    iszero(c) && return T(eigenvalue_Schwarzschild(s, l))
    return N == -1 ? _adaptive_real_lambda(c, s, l, m).value :
        _real_lambda_eigenvalue_at_size(c, s, l, m, N)
end

function Teukolsky_lambda_const(c, s::Int, l::Int, m::Int, N::Int=-1)
    if c isa Complex
        iszero(imag(c)) &&
            return Teukolsky_lambda_const(real(c), s, l, m, N)
        truncation_order = _complex_truncation_order(s, l, m, N, c)
        try
            return continue_angular_mode(s, l, m, c; truncation_order).lambda
        catch exception
            # Where the mode is nearly degenerate with a neighbour, Float64 cannot tell their
            # eigenvectors apart, so the path cannot be followed by overlap; the eigenvalue is still
            # accurate. Only that case is retried, in the wider precision that takes the path over
            # (see _advance_in_precision!), and the eigenvalue is rounded back
            (exception isa AngularContinuationError && _seed_is_final(real(typeof(c))) &&
                exception.spectral_gap < _eigenvector_gap_tol(Float64)) || rethrow()
            return oftype(c, continue_angular_mode(s, l, m, c; truncation_order,
                precision=SWSH_NEAR_DEGENERATE_PRECISION).lambda)
        end
    end
    angular_sep_const(c, s, l, m, N) + c^2 - 2*m*c
end
