using LinearAlgebra
using LinearAlgebra.BLAS: @blasfunc, BlasInt

const SWSH_LAPACK = LinearAlgebra.LAPACK.liblapack
const SWSH_SMALL_C_LIMIT = 1e-6
const SWSH_EIGENVALUE_ATOL = 5e-14
const SWSH_EIGENVALUE_RTOL = 5e-15
const SWSH_MIN_BUFFER = 10
const SWSH_MAX_REFINEMENTS = 20
const SWSH_EIGENVECTOR_OVERLAP_TOL = 5e-13
const SWSH_EIGENVECTOR_RESIDUAL_TOL = 5e-13
# Target size of the neglected spectral coefficients, used to estimate the matrix size
const SWSH_COEFFICIENT_TAIL_TOL = 1e-14
# Number of trailing coefficients checked to decide whether the truncation is adequate
const SWSH_COEFFICIENT_TAIL_WINDOW = 4

#=
The truncation is adequate when the last few coefficients are below
SWSH_COEFFICIENT_TAIL_TOL: the coefficients fall off faster than exponentially,
so the neglected ones beyond N are smaller still.
=#
function _coefficient_tail(coefficients)
    N = length(coefficients)
    return norm(view(coefficients, max(1, N - SWSH_COEFFICIENT_TAIL_WINDOW + 1):N))
end
#=
Roundoff perturbs an eigenvector by an angle of about eps*||M||/gap, where gap is the
distance to the nearest other eigenvalue. Below this relative gap that angle exceeds the
one allowed by SWSH_EIGENVECTOR_OVERLAP_TOL (1 - overlap ≈ angle^2/2), so the eigenvector
cannot be resolved in Float64 even though the eigenvalue still can.
=#
const SWSH_EIGENVECTOR_GAP_TOL = eps(Float64) / sqrt(2 * SWSH_EIGENVECTOR_OVERLAP_TOL)

# TODO: once Leaver's method is merged, point users to it with a wider precision, e.g. big(c)
function _unresolvable_eigenvector_message(s, l, m, c, relative_gap)
    return "The SWSH eigenvector for (s, l, m) = ($s, $l, $m), c = $c cannot be " *
        "resolved in Float64: the mode is nearly degenerate with a neighbouring mode " *
        "(relative gap $relative_gap, below $SWSH_EIGENVECTOR_GAP_TOL), so the " *
        "harmonic would be an arbitrary mix of the two. The eigenvalue is still " *
        "accurate, see spin_weighted_spheroidal_eigenvalue."
end

function _determine_matrix_size_N(s::Int, l::Int, m::Int, c=0)
    #=
    Estimate a suitable value of N for the spectral decomposition

    N covers every spherical harmonic up to the target mode, plus the harmonics
    beyond it that the mode still needs, plus a 'buffer' of SWSH_MIN_BUFFER.
    At large |c| the mode is confined near a pole with a width ~ 1/sqrt(|c|), so its
    spherical-harmonic coefficients fall off like a Gaussian in l of width ~ sqrt(|c|):
    about sqrt(2 log(1/tol) |c|) more harmonics bring them below tol. This is an
    estimate; callers that choose N themselves check the coefficient tail afterwards.
    =#
    N = l - max(abs(m), abs(s)) + 1
    beyond = ceil(Int, sqrt(2log(1 / SWSH_COEFFICIENT_TAIL_TOL) * abs(c)))
    return N + beyond + SWSH_MIN_BUFFER
end

function Fslm(s::Int, l::Int, m::Int)
    # 'Edge' case where l is -1, this can happen when both |m| and |s| are 0 (since lmin = max(|m|, |s|))
    (l == -1 && abs(m) == 0 && abs(s) == 0) && return 0.0
    lp1 = l + 1
    return sqrt((((lp1^2 - m^2) / ((2l + 3) * (2l + 1))) *
        ((lp1^2 - s^2) / lp1^2)))
end

function Gslm(s::Int, l::Int, m::Int)
    iszero(l) && return 0.0
    return sqrt(((l^2 - m^2) / (4l^2 - 1)) *
        ((l^2 - s^2) / l^2))
end

function Hslm(s::Int, l::Int, m::Int)
    (l == 0 || s == 0) ? 0 : -m*s/(l*(l+1))
end

function Aslm(s::Int, l::Int, m::Int)
    Fslm(s, l, m)*Fslm(s, l+1, m)
end

function Bslm(s::Int, l::Int, m::Int)
    Fslm(s, l, m)*Gslm(s, l+1, m) + Gslm(s, l, m)*Fslm(s, l-1, m) + Hslm(s, l, m)^2
end

function Cslm(s::Int, l::Int, m::Int)
    Gslm(s, l, m)*Gslm(s, l-1, m)
end

function Dslm(s::Int, l::Int, m::Int)
    Fslm(s, l, m)*(Hslm(s, l+1, m) + Hslm(s, l, m))
end

function Eslm(s::Int, l::Int, m::Int)
    Gslm(s, l, m)*(Hslm(s, l-1, m) + Hslm(s, l, m))
end

function eigenvalue_Schwarzschild(s::Int, l::Int)
    l*(l+1) - s*(s+1)
end

@inline function _small_c_lambda(c::Real, s::Int, l::Int, m::Int)
    c64 = Float64(c)
    g_lower = Gslm(s, l, m)
    g_upper = Gslm(s, l + 1, m)
    h = Hslm(s, l, m)
    lower_shift = iszero(l) ? 0.0 : g_lower^2 / l
    lambda0 = Float64(eigenvalue_Schwarzschild(s, l))
    lambda1 = 2s * h - 2m
    lambda2 = 1 - Bslm(s, l, m) +
        2s^2 * (lower_shift - g_upper^2 / (l + 1))
    return muladd(c64, muladd(c64, lambda2, lambda1), lambda0)
end

function spectral_matrix_coefficient(c, s::Int, m::Int, l::Int, lprime::Int)
    # This is Eq (55)
    if lprime == l-2
        return -c^2 * Aslm(s, lprime, m)
    elseif lprime == l-1
        return -c^2 * Dslm(s, lprime, m) + 2*c*s*Fslm(s, lprime, m)
    elseif lprime == l
        return eigenvalue_Schwarzschild(s, lprime) - c^2 * Bslm(s, lprime, m) + 2*c*s*Hslm(s, lprime, m)
    elseif lprime == l+1
        return -c^2 * Eslm(s, lprime, m) + 2*c*s*Gslm(s, lprime, m)
    elseif lprime == l+2
        return -c^2 * Cslm(s, lprime, m)
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
    spectral_matrix = zeros(ComplexF64, N, N)
    for i in lmin:lmax
        for j in i:lmax
            spectral_matrix[i-lmin+1, j-lmin+1] = spectral_matrix_coefficient(c, s, m, i, j)
        end
    end
    spectral_matrix = Array(Symmetric(spectral_matrix))
end

function _real_lambda_band(c::Real, s::Int, m::Int, N::Int)
    lmin = max(abs(m), abs(s))
    lmax = lmin + N - 1
    bandwidth = 2
    band = zeros(Float64, bandwidth + 1, N)
    c64 = Float64(c)
    @inbounds for l in lmin:lmax
        i = l - lmin + 1
        for lprime in l:min(l + bandwidth, lmax)
            j = lprime - lmin + 1
            value = if lprime == l
                spherical = Float64(eigenvalue_Schwarzschild(s, l))
                muladd(c64^2, 1 - Bslm(s, l, m),
                    muladd(2c64, s * Hslm(s, l, m) - m, spherical))
            else
                real(spectral_matrix_coefficient(c64, s, m, l, lprime))
            end
            band[bandwidth + 1 + i - j, j] = value
        end
    end
    return band
end

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
function _symmetric_band_inf_norm(band::Matrix{Float64})
    kd = size(band, 1) - 1
    n = size(band, 2)
    row_sums = zeros(Float64, n)
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
Also return the relative gap to the neighbouring eigenvalues. With check_gap, throw if
the eigenvector cannot be resolved (see SWSH_EIGENVECTOR_GAP_TOL).
=#
function _real_lambda_eigenpair_at_size(
    c::Real,
    s::Int,
    l::Int,
    m::Int,
    N::Int;
    check_gap::Bool=true,
)
    index = _ell_index_in_matrix(s, l, m, N)
    lambda, coefficients, relative_gap = _selected_banded_eigenpair!(
        _real_lambda_band(c, s, m, N), index; neighbours=true)
    check_gap && relative_gap < SWSH_EIGENVECTOR_GAP_TOL &&
        error(_unresolvable_eigenvector_message(s, l, m, c, relative_gap))
    return lambda, coefficients, relative_gap
end

function _real_lambda_residual(c::Real, s::Int, m::Int, lambda,
        coefficients)
    N = length(coefficients)
    lmin = max(abs(m), abs(s))
    lmax = lmin + N - 1
    product = zeros(ComplexF64, N)
    maximum_row_sum = 0.0
    c64 = Float64(c)
    @inbounds for l in lmin:lmax
        i = l - lmin + 1
        row_sum = 0.0
        for lprime in l:min(l + 2, lmax)
            j = lprime - lmin + 1
            value = if lprime == l
                spherical = Float64(eigenvalue_Schwarzschild(s, l))
                muladd(c64^2, 1 - Bslm(s, l, m),
                    muladd(2c64, s * Hslm(s, l, m) - m, spherical))
            else
                real(spectral_matrix_coefficient(c64, s, m, l, lprime))
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
    return residual / max(scale, floatmin(Float64))
end

#=
Start from the size estimated from |c|, which is usually enough, and accept the first
solve whose coefficient tail and residual are small: once the neglected coefficients are
below tolerance, a larger matrix could only change the eigenpair at that level. Increase
the size only when the tail is still too large.
=#
function _adaptive_real_eigenpair(c::Real, s::Int, l::Int, m::Int)
    if iszero(c)
        size = _determine_matrix_size_N(s, l, m)
        index = _ell_index_in_matrix(s, l, m, size)
        coefficients = zeros(ComplexF64, size)
        coefficients[index] = 1.0 + 0.0im
        return (lambda=Float64(eigenvalue_Schwarzschild(s, l)),
            coefficients, size, refinement=0, tail=0.0, residual=0.0, gap=Inf)
    end

    c64 = Float64(c)
    size = _determine_matrix_size_N(s, l, m, c64)
    step = max(8, ceil(Int, abs(c64) / 8))
    gap = Inf
    for refinement in 0:SWSH_MAX_REFINEMENTS
        lambda, coefficients, gap = _real_lambda_eigenpair_at_size(
            c64, s, l, m, size; check_gap=false)
        tail = _coefficient_tail(coefficients)
        residual = _real_lambda_residual(c64, s, m, lambda, coefficients)
        if tail <= SWSH_COEFFICIENT_TAIL_TOL && residual <= SWSH_EIGENVECTOR_RESIDUAL_TOL
            gap >= SWSH_EIGENVECTOR_GAP_TOL ||
                error(_unresolvable_eigenvector_message(s, l, m, c, gap))
            return (; lambda, coefficients, size, refinement, tail, residual, gap)
        end
        size += step
    end
    # An unresolvable eigenvector also shows up as one that never settles
    gap >= SWSH_EIGENVECTOR_GAP_TOL ||
        error(_unresolvable_eigenvector_message(s, l, m, c, gap))
    error("SWSH eigenpair failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

function _real_lambda_eigenvalue_at_size(
    c::Real,
    s::Int,
    l::Int,
    m::Int,
    N::Int,
)
    index = _ell_index_in_matrix(s, l, m, N)
    return _selected_banded_eigenvalue!(
        _real_lambda_band(c, s, m, N), index)
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
        return pair.angular_sep, pair.coefficients
    end
    N == -1 || return _dense_real_eigenpair_at_size(c, s, l, m, N)
    N = _determine_matrix_size_N(s, l, m, c)
    step = max(8, ceil(Int, abs(c) / 8))
    for _ in 0:SWSH_MAX_REFINEMENTS
        angular_sep, coefficients = _dense_real_eigenpair_at_size(c, s, l, m, N)
        _coefficient_tail(coefficients) <= SWSH_COEFFICIENT_TAIL_TOL &&
            return angular_sep, coefficients
        N += step
    end
    error("SWSH dense eigenpair did not reach a small enough coefficient tail after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

function _dense_real_eigenpair_at_size(c, s::Int, l::Int, m::Int, N::Int)
    idx = _ell_index_in_matrix(s, l, m, N)
    if c == 0
        coefficients = zeros(ComplexF64, N)
        coefficients[idx] = 1.0 + 0.0im
        return eigenvalue_Schwarzschild(s, l), coefficients
    end
    spectral_matrix = construct_spectral_matrix(c, s, m, N)
    decomposition = eigen(spectral_matrix)
    relative_gap = _relative_spectral_gap(
        decomposition.values, idx, opnorm(spectral_matrix, Inf))
    relative_gap >= SWSH_EIGENVECTOR_GAP_TOL ||
        error(_unresolvable_eigenvector_message(s, l, m, c, relative_gap))
    angular_sep = decomposition.values[idx]
    vector = decomposition.vectors[:, idx]
    pivot = vector[idx]
    !iszero(pivot) && (vector /= pivot)
    coefficients = vector / sqrt(dot(vector, vector))
    return angular_sep, coefficients
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
        angular_sep = lambda -
            muladd(Float64(real(c)), Float64(real(c)),
                -2m * Float64(real(c)))
        return angular_sep, coefficients
    end
    return _dense_spectral_decomposition(c, s, l, m, N)
end

function _adaptive_real_lambda(c::Real, s::Int, l::Int, m::Int)
    iszero(c) && return (value=Float64(eigenvalue_Schwarzschild(s, l)),
        size=0, refinement=0, delta=0.0)
    if abs(c) <= SWSH_SMALL_C_LIMIT
        return (value=_small_c_lambda(c, s, l, m),
            size=_determine_matrix_size_N(s, l, m), refinement=0, delta=0.0)
    end
    lmin = max(abs(m), abs(s))
    idx = l - lmin + 1
    c64 = Float64(c)
    size = max(idx + SWSH_MIN_BUFFER, idx + ceil(Int, abs(c64)) + 8)
    step = max(8, ceil(Int, abs(c64) / 8))
    previous = _real_lambda_eigenvalue_at_size(c64, s, l, m, size)
    for refinement in 1:SWSH_MAX_REFINEMENTS
        next_size = size + step
        current = _real_lambda_eigenvalue_at_size(
            c64, s, l, m, next_size)
        delta = abs(current - previous)
        threshold = SWSH_EIGENVALUE_ATOL +
            SWSH_EIGENVALUE_RTOL * max(abs(previous), abs(current))
        delta <= threshold && return (value=current, size=next_size,
            refinement=refinement, delta=delta)
        size = next_size
        previous = current
    end
    error("SWSH eigenvalue failed to converge after $(SWSH_MAX_REFINEMENTS) matrix refinements.")
end

function _angular_eigenvalue(c::Real, s::Int, l::Int, m::Int, N::Int=-1)
    iszero(c) && return eigenvalue_Schwarzschild(s, l)
    lambda = N == -1 ? _adaptive_real_lambda(c, s, l, m).value :
        _real_lambda_eigenvalue_at_size(c, s, l, m, N)
    return lambda - muladd(Float64(c), Float64(c), -2m * Float64(c))
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
        return continue_angular_mode(
            s, l, m, c;
            truncation_order=_complex_truncation_order(s, l, m, N, c)).coefficients
    end
    _, coeffs = _spectral_decomposition(
        c, s, l, m, N; backend=selected_backend)
    return coeffs
end

function Teukolsky_lambda_const(c::Real, s::Int, l::Int, m::Int, N::Int=-1)
    iszero(c) && return eigenvalue_Schwarzschild(s, l)
    return N == -1 ? _adaptive_real_lambda(c, s, l, m).value :
        _real_lambda_eigenvalue_at_size(c, s, l, m, N)
end

function Teukolsky_lambda_const(c, s::Int, l::Int, m::Int, N::Int=-1)
    if c isa Complex
        iszero(imag(c)) &&
            return Teukolsky_lambda_const(real(c), s, l, m, N)
        return continue_angular_mode(
            s, l, m, c;
            truncation_order=_complex_truncation_order(s, l, m, N, c)).lambda
    end
    angular_sep_const(c, s, l, m, N) + c^2 - 2*m*c
end
