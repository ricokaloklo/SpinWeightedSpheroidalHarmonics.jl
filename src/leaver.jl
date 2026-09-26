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

In contrast to the spectral decomposition, nothing here requires the harmonic
index `l` to be an integer: `l` enters only through the branch index of the
continued fraction and through the initial guess for λ. Note however that the
continued fraction itself keeps a discrete spectrum, so a non-integer `l` does
not by itself continue the eigenvalue; for that, λ has to be supplied (or the
solver started) by hand.
=#

mutable struct LeaverAngularSolution
    c
    lambda
    branch_n::Int
    pminus::Int # |m - s|, so that u^{k₋} = cos(θ/2)^pminus
    pplus::Int  # |m + s|, so that (1 - u)^{k₊} = sin(θ/2)^pplus
    coeffs      # series coefficients a_n, normalized to a_0 = 1
    theta_fun   # Chebyshev representation of the unnormalized S(θ), built on demand for θ-derivatives
end

_LEAVER_TINY = 1e-30 # Continued-fraction denominators smaller than this are nudged away from zero
_LEAVER_CONTINUATION_STEP = 0.25 # Largest step in c taken when marching the eigenvalue away from c = 0
_LEAVER_MAX_UNCERTAINTY = 1e-6 # Relative uncertainty of the eigenvalue beyond which the working precision is deemed insufficient
_LEAVER_MAX_SERIES_ERROR = 1e-3 # Estimated relative accuracy of the power series below which the harmonic is no longer worth returning

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
    return abs(denom) < _LEAVER_TINY ? denom + _LEAVER_TINY * (1 + 1im) : denom
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

    λ = isnothing(lambda0) ? l * (l + 1) - s * (s + 1) : lambda0

    #=
    Newton's method only converges locally, and the eigenvalue drifts by more
    than the spacing between neighboring branches once |c| grows to order unity,
    so a solve started from the c = 0 guess happily converges onto the wrong
    branch. Unless the user supplies their own guess, we therefore march c from
    0 (where the guess above is exact for integer l) to its target value in
    small steps, extrapolating linearly from the two previous eigenvalues.
    =#
    nsteps = isnothing(lambda0) ? max(1, ceil(Int, abs(c) / _LEAVER_CONTINUATION_STEP)) : 1
    λ_previous = λ
    local mismatch
    for k in 1:nsteps
        ck = c * (k / nsteps)
        mismatch = x -> _leaver_mismatch(x, s, m, ck, n; depth=depth)
        λ_solved, converged = _leaver_newton(mismatch, 2 * λ - λ_previous; tol=tolerance, max_iter=max_iter, fd_eps=fd_eps)
        converged || error("Leaver continued-fraction solve did not converge for (s, l, m) = ($s, $l, $m), c = $ck, branch n = $n. Try again with c in a wider floating-point precision, e.g. big($c).")
        λ_previous, λ = λ, λ_solved
    end

    uncertainty = _leaver_root_uncertainty(mismatch, λ, working_eps)
    uncertainty > _LEAVER_MAX_UNCERTAINTY * (1 + abs(λ)) && error("Roundoff in the continued fraction only locates the eigenvalue for (s, l, m) = ($s, $l, $m), c = $c, branch n = $n to within $uncertainty. Try again with c in a wider floating-point precision, e.g. big($c).")
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
    _leaver_converged_series_coefficients(s::Int, m::Int, c, λ, branch_n::Int; nmax::Int=-1)

Return the series coefficients `a_n` truncated where they stop contributing at
machine precision. When `nmax = -1`, the number of terms is grown automatically
until the tail has decayed; otherwise exactly `nmax + 1` terms are computed.
"""
function _leaver_converged_series_coefficients(s::Int, m::Int, c, λ, branch_n::Int; nmax::Int=-1, nmax_limit::Int=6400)
    # The branch index has to be inside the series for the two halves of the
    # recurrence in _leaver_series_coefficients to meet
    smallest_nmax = branch_n + 1
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

@doc raw"""
    _leaver_raw_theta_value(sol::LeaverAngularSolution, theta)

Evaluate the *unnormalized* θ-dependent part of the harmonic. `theta` must
already lie in ``[0, \pi]``.

We write ``u^{k_-} (1-u)^{k_+} = \cos(\theta/2)^{|m-s|} \sin(\theta/2)^{|m+s|}``,
which keeps the exponents integral and makes the manifest smoothness of `S` in
θ (as opposed to in ``\cos\theta``) available to the Chebyshev representation
used for θ-derivatives.
"""
function _leaver_prefactor(sol::LeaverAngularSolution, theta)
    half_theta = theta / 2
    cos_half, sin_half = cos(half_theta), sin(half_theta)
    u = cos_half^2
    return exp(2 * sol.c * u) * cos_half^sol.pminus * sin_half^sol.pplus, u
end

function _leaver_raw_theta_value(sol::LeaverAngularSolution, theta)
    prefactor, u = _leaver_prefactor(sol, theta)
    return prefactor * _leaver_series_sum(sol.coeffs, u)
end

#=
Fit `f` on [0, π] with a Chebyshev series of fixed length, doubling the length
until its tail has decayed. Evaluating the Leaver series involves cancellations
that leave a roundoff floor, and letting ApproxFun construct the Fun adaptively
would make it refine against that floor forever instead of stopping.
=#
function _leaver_chebyshev_fit(f; n_start::Int=64, n_max::Int=4096)
    space = Chebyshev(0..π)
    n = n_start
    while true
        coefficients = transform(space, f.(points(space, n)))
        rtol = 100 * eps(real(eltype(coefficients)))
        tail = maximum(abs, @view coefficients[end-n÷8:end])
        (tail <= rtol * maximum(abs, coefficients) || n >= n_max) && return Fun(space, coefficients)
        n *= 2
    end
end

function _leaver_theta_fun(sol::LeaverAngularSolution)
    if isnothing(sol.theta_fun)
        sol.theta_fun = _leaver_chebyshev_fit(theta -> _leaver_raw_theta_value(sol, theta))
    end
    return sol.theta_fun
end

@doc raw"""
    _leaver_normalization_const(sol::LeaverAngularSolution, s::Int, l, m::Int)

Compute the complex constant `K` such that `S = S_raw / K` follows the same
convention as the spectral decomposition, namely
``\int_0^\pi |S(\theta, \phi)|^2 \sin\theta \, d\theta = 1/(2\pi)`` with the
overlap ``\langle {}_s Y_{lm} | S \rangle`` real and positive.

The phase can only be fixed this way when `l` is an integer; for a continued
(non-integer) `l` there is no spherical harmonic to project onto and we instead
keep the natural Leaver convention `a_0 = 1`.
"""
function _leaver_normalization_const(sol::LeaverAngularSolution, s::Int, l, m::Int)
    theta_integral = sum(_leaver_chebyshev_fit(theta -> abs2(_leaver_raw_theta_value(sol, theta)) * sin(theta)))
    theta_integral <= 0 && error("Could not normalize the Leaver spheroidal harmonic: the norm is non-positive.")
    magnitude = sqrt(2π * theta_integral)

    if l isa Integer && l >= max(abs(s), abs(m))
        spherical = spin_weighted_spherical_harmonic(s, Int(l), m)
        overlap = sum(_leaver_chebyshev_fit(theta -> conj(spherical(theta, 0.0)) * _leaver_raw_theta_value(sol, theta) * sin(theta)))
        # |overlap| is bounded by sqrt(theta_integral/(2π)); bail out of fixing
        # the phase if this mode happens to have essentially no overlap with sYlm
        if abs(overlap) > 1e-3 * sqrt(theta_integral / (2π))
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
    sol = LeaverAngularSolution(c, λ, branch_n, abs(m - s), abs(m + s), coeffs, nothing)

    accuracy = _leaver_series_accuracy(sol)
    accuracy > _LEAVER_MAX_SERIES_ERROR && error("Cancellation in Leaver's power series for (s, m) = ($s, $m), c = $c, branch n = $branch_n leaves it with a relative accuracy of only about $accuracy. Use the spectral decomposition (e.g. method=\"jacobi\"), which does not suffer from this, or a wider floating-point precision.")

    return sol
end

function _nth_derivative_spheroidal_harmonic_leaver(sol::LeaverAngularSolution, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    _theta, _phi = _fold_theta_phi(theta, phi)

    # Note that the derivatives at the two boundary points might be inaccurate
    S_deriv = theta_derivative == 0 ? _leaver_raw_theta_value(sol, _theta) : differentiate(_leaver_theta_fun(sol), theta_derivative)(_theta)
    return S_deriv * cis(m * _phi) * (m * 1im)^phi_derivative
end
