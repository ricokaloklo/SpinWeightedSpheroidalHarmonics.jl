#=
Leaver's continued-fraction method for the angular Teukolsky equation
(E. W. Leaver, Proc. R. Soc. Lond. A 402, 285 (1985)).

The angular function is written as

    S(θ) = e^{2cu} u^{k₋} (1 - u)^{k₊} Σ_{n ≥ 0} a_n u^n ,   u = (1 + cosθ)/2 ,

with k₋ = |m - s|/2 and k₊ = |m + s|/2. The coefficients obey the three-term
recurrence

    α_n a_{n+1} + β_n a_n + γ_n a_{n-1} = 0 ,

which admits a minimal (i.e. normalizable) solution only for the eigenvalues λ.
Those are therefore the roots of a continued fraction, which we invert with a
damped Newton iteration.

In the eigenvalue solve, the harmonic index `l` enters only through the branch index of the
continued fraction and through the initial guess for λ. The public API requires an integer
`l` (see `_resolve_spheroidal_method`); a non-integer degree needs a dedicated solver.
=#

struct LeaverAngularSolution
    c
    lambda
    branch_n::Int
    pminus::Int # |m - s|, so that u^{k₋} = cos(θ/2)^pminus
    pplus::Int  # |m + s|, so that (1 - u)^{k₊} = sin(θ/2)^pplus
    coeffs      # series coefficients a_n, normalized to a_0 = 1
end

# Continued-fraction denominators smaller than this are nudged away from zero
_leaver_tiny(::Type{T}) where {T} = _tol(T, 1e-30)
const _LEAVER_CONTINUATION_STEP = 0.25 # Largest step in c taken when marching the eigenvalue away from c = 0
# Relative uncertainty of the eigenvalue beyond which the working precision is deemed insufficient
_leaver_max_uncertainty(::Type{T}) where {T} = _tol(T, 1e-6)
# Estimated relative error of the power series beyond which the harmonic is no longer worth returning
_leaver_max_series_error(::Type{T}) where {T} = _tol(T, 1e-3)

@doc raw"""
    _leaver_recurrence_coefficients(n::Int, s::Int, m::Int, c, λ)

Return the three-term recurrence coefficients `(α_n, β_n, γ_n)` of Leaver's
continued-fraction method for the angular Teukolsky equation, i.e.
$\alpha_n a_{n+1} + \beta_n a_n + \gamma_n a_{n-1} = 0$.
"""
function _leaver_recurrence_coefficients(n::Int, s::Int, m::Int, c, λ)
    kminus = abs(m - s) / 2
    kplus = abs(m + s) / 2

    αn = (n + 1) * (n + 2 * kminus + 1)
    β0 = λ + s * (s + 1) - m * s +
         2 * c * (2 * kminus + m + s + 1) -
         (2 * kminus^2 + 2 * kminus * kplus + kminus + kplus)
    βn = -n * (n - 1) + (4 * c - 2 * kminus - 2 * kplus - 2) * n + β0
    γn = -4 * c * (n + kminus + kplus + s)
    return αn, βn, γn
end

function _leaver_branch_index(s::Int, m::Int, l)
    lmin = max(abs(s), abs(m))
    n = round(real(l - lmin))
    if !isfinite(n) || n < 0
        error("Cannot infer a Leaver branch index from l=$l (with s=$s, m=$m, so lmin=$lmin). Pass branch_n explicitly.")
    end
    return Int(n)
end

function _leaver_regularize_denominator(denom)
    tiny = _leaver_tiny(real(typeof(float(denom))))
    return abs(denom) < tiny ? denom + tiny * (1 + 1im) : denom
end

# Continued fraction running from the n-th coefficient down to a_0
function _leaver_left_cf_term(s::Int, m::Int, c, λ, n::Int)
    n == 0 && return zero(λ + c)

    _, β0, _ = _leaver_recurrence_coefficients(0, s, m, c, λ)
    denom = β0
    for k in 1:n-1
        αkm1, _, _ = _leaver_recurrence_coefficients(k - 1, s, m, c, λ)
        _, βk, γk = _leaver_recurrence_coefficients(k, s, m, c, λ)
        denom = βk - αkm1 * γk / _leaver_regularize_denominator(denom)
    end

    αnm1, _, _ = _leaver_recurrence_coefficients(n - 1, s, m, c, λ)
    _, _, γn = _leaver_recurrence_coefficients(n, s, m, c, λ)
    return αnm1 * γn / _leaver_regularize_denominator(denom)
end

# Continued fraction running from the n-th coefficient up to the truncated tail
function _leaver_right_cf_term(s::Int, m::Int, c, λ, n::Int; depth::Int=1000)
    nmax = n + depth
    _, βtail, _ = _leaver_recurrence_coefficients(nmax, s, m, c, λ)
    denom = βtail

    for k in (nmax-1):-1:(n+1)
        αk, βk, _ = _leaver_recurrence_coefficients(k, s, m, c, λ)
        _, _, γkp1 = _leaver_recurrence_coefficients(k + 1, s, m, c, λ)
        denom = βk - αk * γkp1 / _leaver_regularize_denominator(denom)
    end

    αn, _, _ = _leaver_recurrence_coefficients(n, s, m, c, λ)
    _, _, γnp1 = _leaver_recurrence_coefficients(n + 1, s, m, c, λ)
    return αn * γnp1 / _leaver_regularize_denominator(denom)
end

@doc raw"""
    _leaver_mismatch(λ, s::Int, m::Int, c, n::Int; depth::Int=1000)

Compute the residual of the `n`-th inversion of the continued fraction, i.e. the
condition that the two halves of the recurrence agree at `a_n`. Its roots are
the eigenvalues λ, and they are the same for every `n`; what changes with `n` is
how well conditioned the root-finding is.
"""
function _leaver_mismatch(λ, s::Int, m::Int, c, n::Int; depth::Int=1000)
    _, βn, _ = _leaver_recurrence_coefficients(n, s, m, c, λ)
    return βn - _leaver_left_cf_term(s, m, c, λ, n) - _leaver_right_cf_term(s, m, c, λ, n; depth=depth)
end

# Precision of the arithmetic the continued fraction will be evaluated in
_leaver_working_eps(c) = eps(typeof(real(float(one(c)))))

# Damped Newton iteration with a central finite-difference derivative.
# The residual is analytic in λ but has poles in between its roots, hence the damping.
function _leaver_newton(f, λ0; tol=1e-12, max_iter::Int=80, fd_eps=1e-7)
    λ = λ0
    best_λ, best_residual = λ0, Inf
    stalled = 0

    for _ in 1:max_iter
        F = f(λ)
        abs(F) <= tol * (1 + abs(λ)) && return λ, true

        if abs(F) < best_residual
            best_λ, best_residual = λ, abs(F)
            stalled = 0
        else
            stalled += 1
        end

        h = fd_eps * (1 + abs(λ))
        dF = (f(λ + h) - f(λ - h)) / (2 * h)
        if abs(dF) == 0
            # Probe along an independent direction in the complex plane
            h *= (1 + 1im)
            dF = (f(λ + h) - f(λ - h)) / (2 * h)
            abs(dF) == 0 && return best_λ, false
        end

        step = F / dF
        #=
        Cancellation between the continued fraction and β_n puts a roundoff
        floor under the residual, which for large |c| sits well above `tol`.
        Once the iterate has settled and the residual has stopped improving,
        that floor is all that is left, so take the best iterate seen.
        =#
        if stalled >= 3 && abs(step) <= sqrt(tol) * (1 + abs(λ))
            return best_λ, true
        end

        trial = λ - step
        Ftrial = f(trial)
        # Halve the step at most six times; after that the trial point is taken even if the residual grew
        damping = 1.0
        while abs(Ftrial) > abs(F) && damping > 1 / 64
            damping /= 2
            trial = λ - damping * step
            Ftrial = f(trial)
        end
        λ = trial
    end
    return best_λ, false
end

#=
Estimate how sharply the working precision pins down the root. Cancellation
between the continued fraction and β_n leaves a roundoff floor on the residual,
and the root is only located to within that floor divided by the local slope.
The floor grows steeply with |c|: in Float64 it costs nothing up to |c| ~ 4,
but by |c| ~ 12 it leaves barely two significant digits of the eigenvalue.
=#
function _leaver_root_uncertainty(f, λ, working_eps)
    scale = 1 + abs(λ)
    # The spread of the residual over seven points a thousand ulps apart around λ is taken as its roundoff floor
    jitter = 1000 * working_eps * scale
    residuals = [abs(f(λ + k * jitter)) for k in -3:3]
    noise = maximum(residuals) - minimum(residuals)

    h = cbrt(working_eps) * scale
    slope = abs(f(λ + h) - f(λ - h)) / (2 * h)
    return slope == 0 ? typeof(noise)(Inf) : noise / slope
end

# The continued fraction settles within a few tens of levels over the range of
# spheroidicities where a Float64 evaluation of it means anything at all, so
# this keeps a comfortable margin on top of what is actually needed
_leaver_default_cf_depth(n::Int, c) = max(200, 20 * n + 50 * ceil(Int, abs(c)))

const _LEAVER_MIN_OVERLAP = 0.9 # Smallest overlap between the eigenfunctions of consecutive march steps for a step to be accepted
const _LEAVER_MAX_HALVINGS = 12 # A march step may shrink to 2^-12 of the nominal step before the march gives up
const _LEAVER_MARCH_BUDGET = 10 # The march gives up after this many solves per nominal step (plus a fixed 100)
# θ grid on which _leaver_mode_profile samples the eigenfunction; the poles, where the sqrt(sin θ) weight
# vanishes, are left out
const _LEAVER_PROFILE_THETA = range(0.02, π - 0.02, length=48)

# `where` is the c at which the solve failed, which during the march in c can differ from the target c
_leaver_nonconvergence_message(s, l, m, c, n; where=c) = "Leaver continued-fraction solve did not converge for (s, l, m) = ($s, $l, $m), c = $where, branch n = $n. Try again with c in a wider floating-point precision, e.g. big($c)."

# The eigenfunction for λ, sampled on a θ grid (weighted by sqrt(sin θ)) and scaled to unit norm,
# so that the modes found at two nearby c can be compared
function _leaver_mode_profile(s::Int, m::Int, c, λ, n::Int)
    coeffs = _leaver_converged_series_coefficients(s, m, c, λ, n)
    sol = LeaverAngularSolution(c, λ, n, abs(m - s), abs(m + s), coeffs)
    values = [_leaver_raw_theta_value(sol, θ) * sqrt(sin(θ)) for θ in _LEAVER_PROFILE_THETA]
    return values ./ sqrt(sum(abs2, values))
end

#=
March c from 0 to its target along a straight line, starting from the eigenvalue λ0 at c = 0.

Following λ alone is not enough: once a neighbouring mode comes close, the extrapolated
guess can lie nearer that mode's root, and Newton converges there, in any precision. So each
step is also checked against the eigenfunction. Between nearby c the same mode changes very
little, while a neighbouring mode is a different function (for real c, an orthogonal one).
A step that fails either check is halved, and steps grow back once the march is past the
difficulty. Where the neighbour is so close that only tiny steps would stay on the mode,
the march gives up rather than grinding on: both the smallest step and the total number of
solves are limited. Return the eigenvalue and the mismatch function (the residual in λ) at the target c.
=#
function _leaver_march(s::Int, l, m::Int, c, n::Int, λ0; depth, tolerance, max_iter, fd_eps)
    nominal = iszero(c) ? 1.0 : min(1.0, _LEAVER_CONTINUATION_STEP / Float64(abs(c)))
    smallest_step = nominal / 2^_LEAVER_MAX_HALVINGS
    budget = _LEAVER_MARCH_BUDGET * ceil(Int, 1 / nominal) + 100
    t_previous, λ_previous = 0.0, λ0
    t, λ = 0.0, λ0
    profile = _leaver_mode_profile(s, m, zero(c), λ0, n)
    step = nominal
    solves = 0
    failure = :none
    mismatch = x -> _leaver_mismatch(x, s, m, zero(c), n; depth=depth)
    give_up() = failure == :newton ? error(_leaver_nonconvergence_message(s, l, m, c, n; where=c * t)) :
        error("Leaver's continued fraction cannot follow the mode (s, l, m) = ($s, $l, $m) " *
            "past c = $(c * t): it becomes nearly degenerate with a neighbouring mode, and the " *
            "solve keeps converging onto that mode. Seed it with the spectral eigenvalue instead, " *
            "lambda0 = spin_weighted_spheroidal_eigenvalue($s, $l, $m, $(isreal(c) ? Float64(real(c)) : ComplexF64(c))).")
    while t < 1
        solves += 1
        solves > budget && give_up()
        t_next = min(1.0, t + step)
        c_next = c * t_next
        # Linear extrapolation from the last two accepted points (the steps can differ in size)
        guess = t == t_previous ? λ : λ + (λ - λ_previous) * ((t_next - t) / (t - t_previous))
        f = x -> _leaver_mismatch(x, s, m, c_next, n; depth=depth)
        λ_next, converged = _leaver_newton(f, guess; tol=tolerance, max_iter=max_iter, fd_eps=fd_eps)
        next_profile = nothing
        if !converged
            failure = :newton
        else
            next_profile = try
                _leaver_mode_profile(s, m, c_next, λ_next, n)
            catch
                nothing # λ_next is not an eigenvalue that gives a normalizable solution
            end
            if next_profile === nothing || abs(sum(conj.(profile) .* next_profile)) < _LEAVER_MIN_OVERLAP
                failure = :mode
                next_profile = nothing
            end
        end
        if next_profile !== nothing
            t_previous, λ_previous = t, λ
            t, λ, profile, mismatch = t_next, λ_next, next_profile, f
            step = min(2step, nominal)
        else
            step /= 2
            step < smallest_step && give_up()
        end
    end
    return λ, mismatch
end

@doc raw"""
    _leaver_eigenvalue(s::Int, l, m::Int, c; branch_n=nothing, lambda0=nothing, cf_depth::Int=-1, tol=-1, max_iter::Int=80)

Compute the spin-weighted spheroidal eigenvalue λ with Leaver's continued fraction.
"""
function _leaver_eigenvalue(s::Int, l, m::Int, c;
    branch_n::Union{Nothing,Int}=nothing,
    lambda0=nothing,
    cf_depth::Int=-1,
    tol=-1,
    max_iter::Int=80)

    n = isnothing(branch_n) ? _leaver_branch_index(s, m, l) : branch_n
    n < 0 && error("Branch index n=$n is invalid. Choose branch_n >= 0.")
    depth = cf_depth == -1 ? _leaver_default_cf_depth(n, c) : cf_depth

    # Both the tolerance and the step of the finite-difference derivative follow
    # the precision `c` is given in, so that a wider precision is actually put to
    # use rather than being thrown away by a Float64-sized stopping criterion
    working_eps = _leaver_working_eps(c)
    tolerance = tol == -1 ? 1000 * working_eps : tol
    fd_eps = cbrt(working_eps)

    #=
    Newton's method only converges locally, and the eigenvalue drifts by more
    than the spacing between neighbouring branches once |c| grows to order unity,
    so a solve started from the c = 0 guess happily converges onto the wrong
    branch. Unless the user supplies their own guess, we therefore march c from
    0 (where l(l+1) - s(s+1) is exact for integer l) to its target value.
    =#
    T = real(typeof(float(c)))
    if isnothing(lambda0)
        λ, mismatch = _leaver_march(s, l, m, c, n, T(l * (l + 1) - s * (s + 1));
            depth=depth, tolerance=tolerance, max_iter=max_iter, fd_eps=fd_eps)
    else
        mismatch = x -> _leaver_mismatch(x, s, m, c, n; depth=depth)
        λ, converged = _leaver_newton(mismatch, lambda0; tol=tolerance, max_iter=max_iter, fd_eps=fd_eps)
        converged || error(_leaver_nonconvergence_message(s, l, m, c, n))
    end

    uncertainty = _leaver_root_uncertainty(mismatch, λ, working_eps)
    uncertainty > _leaver_max_uncertainty(T) * (1 + abs(λ)) && error("Roundoff in the continued fraction only locates the eigenvalue for (s, l, m) = ($s, $l, $m), c = $c, branch n = $n to within $uncertainty. Try again with c in a wider floating-point precision, e.g. big($c).")
    return λ
end

#=
Solve the recurrence for the minimal (i.e. normalizable) solution a_n.

The tail comes from a backward sweep for the ratios a_{n+1}/a_n, seeded far
beyond the coefficients we actually keep; recursing forward there would be
swamped by the dominant solution. The backward sweep is however ill-conditioned
below the branch index (for the n-th branch it only reproduces a_0, ..., a_n
through a cancellation that becomes exactly 0/0 as c → 0), so those first
coefficients are instead recursed forward, where the dominant solution has not
had room to grow. That the two halves agree where they meet is precisely the
n-th inversion of the continued fraction, i.e. the condition that fixes λ.
=#
function _leaver_series_coefficients(s::Int, m::Int, c, λ, branch_n::Int, nmax::Int, ratio_tail::Int)
    nmax < 0 && error("nmax must be non-negative.")
    ratio_tail < nmax + 4 && error("ratio_tail must be at least nmax+4 for a stable backward recurrence.")

    T = typeof(complex(λ + c))
    ratios = Vector{T}(undef, nmax + 1) # ratios[n+1] = a_{n+1}/a_n

    tail_ratio = zero(T)
    for k in ratio_tail:-1:0
        αkp1, βkp1, γkp1 = _leaver_recurrence_coefficients(k + 1, s, m, c, λ)
        rk = -γkp1 / _leaver_regularize_denominator(βkp1 + αkp1 * tail_ratio)
        k <= nmax && (ratios[k+1] = rk)
        tail_ratio = rk
    end

    coeffs = Vector{T}(undef, nmax + 1)
    coeffs[1] = one(T)
    for n in 0:nmax-1
        if n < branch_n
            αn, βn, γn = _leaver_recurrence_coefficients(n, s, m, c, λ)
            previous = n == 0 ? zero(T) : coeffs[n]
            coeffs[n+2] = -(βn * coeffs[n+1] + γn * previous) / αn
        else
            coeffs[n+2] = ratios[n+1] * coeffs[n+1]
        end
    end
    return coeffs
end

# Drop the tail of the series that no longer contributes at machine precision
function _leaver_trim_series(coeffs)
    largest = maximum(abs, coeffs)
    largest == 0 && return coeffs
    cutoff = eps(typeof(largest)) * largest
    last_significant = findlast(coeff -> abs(coeff) > cutoff, coeffs)
    return coeffs[1:max(something(last_significant, 1), 2)]
end

@doc raw"""
    _leaver_converged_series_coefficients(s::Int, m::Int, c, λ, branch_n::Int; nmax::Int=-1, nmax_limit)

Return the series coefficients `a_n` truncated where they stop contributing at the working
precision. When `nmax = -1`, the number of terms is grown automatically, up to `nmax_limit`
(6400 in Float64, more in a wider precision), until the tail has decayed; otherwise
`max(nmax, branch_n + 1) + 1` terms are computed and trimmed in the same way.
"""
function _leaver_converged_series_coefficients(s::Int, m::Int, c, λ, branch_n::Int; nmax::Int=-1,
        nmax_limit::Int=_nterms(real(typeof(float(c))), 6400))
    # The branch index has to be inside the series for the two halves of the
    # recurrence in _leaver_series_coefficients to meet
    smallest_nmax = branch_n + 1
    # The backward sweep for the ratios is seeded 800 levels beyond the last computed coefficient, far
    # enough for it to have settled onto the minimal solution by the time it reaches the coefficients
    # that are kept (see _leaver_series_coefficients)
    if nmax != -1
        n = max(nmax, smallest_nmax)
        return _leaver_trim_series(_leaver_series_coefficients(s, m, c, λ, branch_n, n, n + 800))
    end

    n = max(200, 2 * smallest_nmax)
    while n <= nmax_limit
        coeffs = _leaver_series_coefficients(s, m, c, λ, branch_n, n, n + 800)
        trimmed = _leaver_trim_series(coeffs)
        # A shorter result means the tail decayed below machine precision on its own
        length(trimmed) < length(coeffs) && return trimmed
        n *= 2
    end
    error("The Leaver series for (s, m) = ($s, $m), c = $c, branch n = $branch_n did not decay within $nmax_limit terms. This usually means that lambda = $λ is not an eigenvalue of that branch.")
end

function _leaver_series_sum(coeffs, u)
    y = coeffs[end]
    for k in (length(coeffs)-1):-1:1
        y = coeffs[k] + u * y
    end
    return y
end

# The real floating-point type the solution is computed in
_leaver_real_type(sol::LeaverAngularSolution) = real(eltype(sol.coeffs))

function _leaver_prefactor(sol::LeaverAngularSolution, theta)
    half_theta = theta / 2
    cos_half, sin_half = cos(half_theta), sin(half_theta)
    u = cos_half^2
    return exp(2 * sol.c * u) * cos_half^sol.pminus * sin_half^sol.pplus, u
end

@doc raw"""
    _leaver_raw_theta_value(sol::LeaverAngularSolution, theta)

Evaluate the *unnormalized* θ-dependent part of the harmonic. `theta` must
already lie in ``[0, \pi]``.

The prefactor, computed by `_leaver_prefactor`, is written as
``e^{2cu} u^{k_-} (1-u)^{k_+} = e^{2cu} \cos(\theta/2)^{|m-s|} \sin(\theta/2)^{|m+s|}``,
which keeps the exponents integral, so that the harmonic and its θ-derivatives
(see `_leaver_theta_derivative`) are smooth functions of θ, also at the poles.
"""
function _leaver_raw_theta_value(sol::LeaverAngularSolution, theta)
    prefactor, u = _leaver_prefactor(sol, theta)
    return prefactor * _leaver_series_sum(sol.coeffs, u)
end

#=
θ-derivatives of the unnormalized S(θ) = exp(2cu) cos(θ/2)^pminus sin(θ/2)^pplus Σ a_n u^n, u = cos²(θ/2),
from its Taylor coefficients in t at θ = θ0 + t, truncated after t^k ("jets"). Products, powers, the
exponential and the series are exact operations on the coefficients, so the derivatives are as accurate as
the value itself and need no division: they hold at the poles too, where cos(θ/2) or sin(θ/2) vanishes.
The k-th derivative is k! times the k-th coefficient.
=#
_jet_mul(a, b) = [sum(a[i+1] * b[k-i+1] for i in 0:k) for k in 0:length(a)-1]

function _jet_pow(a, n::Int)
    result = [k == 0 ? one(eltype(a)) : zero(eltype(a)) for k in 0:length(a)-1]
    power = a
    while n > 0
        isodd(n) && (result = _jet_mul(result, power))
        n >>= 1
        n > 0 && (power = _jet_mul(power, power))
    end
    return result
end

function _jet_exp(g)
    e = similar(g)
    e[1] = exp(g[1])
    for k in 1:length(g)-1
        e[k+1] = sum(j * g[j+1] * e[k-j+1] for j in 1:k) / k
    end
    return e
end

function _leaver_theta_derivative(sol::LeaverAngularSolution, theta, k::Int)
    T = typeof(float(theta))
    cos_half, sin_half = cos(theta / 2), sin(theta / 2)
    # cos(θ/2 + t/2) and sin(θ/2 + t/2): the j-th derivative of cos is cos, -sin, -cos, sin, … and each brings 1/2
    scale(j) = inv(T(2)^j * prod(T, 1:j; init=one(T)))
    C = [(cos_half, -sin_half, -cos_half, sin_half)[j % 4 + 1] * scale(j) for j in 0:k]
    S = [(sin_half, cos_half, -sin_half, -cos_half)[j % 4 + 1] * scale(j) for j in 0:k]
    u = _jet_mul(C, C)
    series = fill(zero(eltype(sol.coeffs)) * zero(T), k + 1)
    for coefficient in Iterators.reverse(sol.coeffs)
        series = _jet_mul(series, u)
        series[1] += coefficient
    end
    jet = _jet_mul(_jet_mul(_jet_exp(2 * sol.c .* u), _jet_pow(C, sol.pminus)), _jet_mul(_jet_pow(S, sol.pplus), series))
    return prod(T, 1:k; init=one(T)) * jet[k+1]
end

#=
Fit `f` on [0, π] with a Chebyshev series of fixed length, doubling the length
until its tail has decayed. Evaluating the Leaver series involves cancellations
that leave a roundoff floor, and letting ApproxFun construct the Fun adaptively
would make it refine against that floor forever instead of stopping. So the doubling
also stops once the tail no longer falls (by at least a factor of 4) and lies within a
factor of 100 of `floor`, the estimated roundoff of the values of `f` relative to the
largest one (see _leaver_series_accuracy): it has reached that floor, and the longer
of the two fits is kept. (A tail that has not started to fall yet, as for a harmonic of
large index at small n, is far above the floor.)
=#
function _leaver_chebyshev_fit(f, ::Type{T}, floor; n_start::Int=64, n_max::Int=4096) where {T}
    space = Chebyshev(zero(T)..T(π))
    n = n_start
    previous_tail = T(Inf)
    while true
        coefficients = transform(space, f.(points(space, n)))
        rtol = 100 * eps(real(eltype(coefficients)))
        tail = maximum(abs, @view coefficients[end-n÷8:end]) / maximum(abs, coefficients)
        (tail <= rtol || (tail <= 100floor && tail > previous_tail / 4) || n >= n_max) &&
            return Fun(space, coefficients)
        previous_tail = tail
        n *= 2
    end
end

@doc raw"""
    _leaver_normalization_const(sol::LeaverAngularSolution, s::Int, l, m::Int; analytic=false)

Compute the complex constant `K` such that `S = S_raw / K` follows the same
convention as the spectral decomposition, namely
``\int_0^\pi |S(\theta, \phi)|^2 \sin\theta \, d\theta = 1/(2\pi)`` with the
overlap ``\langle {}_s Y_{lm} | S \rangle`` real and positive.

The phase can only be fixed this way when `l` is an integer; for a continued
(non-integer) `l` there is no spherical harmonic to project onto and we instead
keep the natural Leaver convention `a_0 = 1`.

With `analytic=true` and complex `c`, `K` instead makes
``\int_0^\pi S(\theta, 0)^2 \sin\theta \, d\theta = 1/(2\pi)``, the square root taken
so that the overlap with ``{}_s Y_{lm}`` lies in the right half-plane (for real `c` the two agree).
"""
function _leaver_normalization_const(sol::LeaverAngularSolution, s::Int, l, m::Int; analytic::Bool=false)
    T = _leaver_real_type(sol)
    floor = _leaver_series_accuracy(sol)
    if analytic && !isreal(sol.c) && l isa Integer
        squared_integral = sum(_leaver_chebyshev_fit(theta -> _leaver_raw_theta_value(sol, theta)^2 * sin(theta), T, floor))
        scale = sqrt(2T(π) * squared_integral)
        # sYlm(theta, 0) is real, so the overlap needs no conjugation either
        spherical = spin_weighted_spherical_harmonic(s, Int(l), m; precision=_precision_keyword(_precision_bits(T)))
        overlap = sum(_leaver_chebyshev_fit(theta -> real(spherical(theta, zero(T))) * _leaver_raw_theta_value(sol, theta) * sin(theta), T, floor))
        return real(overlap / scale) < 0 ? -scale : scale
    end
    theta_integral = sum(_leaver_chebyshev_fit(theta -> abs2(_leaver_raw_theta_value(sol, theta)) * sin(theta), T, floor))
    theta_integral <= 0 && error("Could not normalize the Leaver spheroidal harmonic: the norm is non-positive.")
    magnitude = sqrt(2T(π) * theta_integral)

    if l isa Integer && l >= max(abs(s), abs(m))
        spherical = spin_weighted_spherical_harmonic(s, Int(l), m; precision=_precision_keyword(_precision_bits(T)))
        overlap = sum(_leaver_chebyshev_fit(theta -> conj(spherical(theta, zero(T))) * _leaver_raw_theta_value(sol, theta) * sin(theta), T, floor))
        # |overlap| is bounded by sqrt(theta_integral/(2π)); bail out of fixing
        # the phase if this mode happens to have essentially no overlap with sYlm
        if abs(overlap) > 1e-3 * sqrt(theta_integral / (2T(π)))
            return magnitude * overlap / abs(overlap)
        end
    end
    return magnitude
end

#=
Estimate the accuracy the series can still deliver, relative to the overall size
of the harmonic, from the largest roundoff its own summation can produce. Powers
of u are a badly conditioned basis for a harmonic of large index: the
coefficients reach ~10^(2l/3) while their sum stays of order unity, and that
cancellation eats into the working precision long before the eigenvalue itself
becomes hard to find.
=#
function _leaver_series_accuracy(sol::LeaverAngularSolution)
    magnitudes = abs.(sol.coeffs)
    largest_error = zero(eltype(magnitudes))
    largest_value = zero(eltype(magnitudes))

    for theta in range(0, π; length=65)
        prefactor, u = _leaver_prefactor(sol, theta)
        # Measuring against the largest value, not against the value here, so
        # that a node of the harmonic does not look like a loss of precision
        largest_error = max(largest_error, abs(prefactor) * _leaver_series_sum(magnitudes, u))
        largest_value = max(largest_value, abs(prefactor * _leaver_series_sum(sol.coeffs, u)))
    end

    largest_value == 0 && return typeof(largest_error)(Inf)
    return largest_error * eps(eltype(magnitudes)) / largest_value
end

@doc raw"""
    _leaver_angular_solution(s::Int, m::Int, c, λ; branch_n::Int, nmax::Int=-1)

Build the representation of the Leaver series for the eigenvalue `λ`.
"""
function _leaver_angular_solution(s::Int, m::Int, c, λ; branch_n::Int, nmax::Int=-1)
    coeffs = _leaver_converged_series_coefficients(s, m, c, λ, branch_n; nmax=nmax)
    sol = LeaverAngularSolution(c, λ, branch_n, abs(m - s), abs(m + s), coeffs)

    accuracy = _leaver_series_accuracy(sol)
    accuracy > _leaver_max_series_error(_leaver_real_type(sol)) && error("Cancellation in Leaver's power series for (s, m) = ($s, $m), c = $c, branch n = $branch_n leaves it with a relative accuracy of only about $accuracy. Use the spectral decomposition (e.g. method=\"jacobi\"), which does not suffer from this, or a wider floating-point precision.")

    return sol
end

function _nth_derivative_spheroidal_harmonic_leaver(sol::LeaverAngularSolution, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    _theta, _phi = _fold_theta_phi(theta, phi)

    S_deriv = theta_derivative == 0 ? _leaver_raw_theta_value(sol, _theta) :
        _leaver_theta_derivative(sol, _theta, theta_derivative)
    return S_deriv * cis(m * _phi) * (m * 1im)^phi_derivative
end
