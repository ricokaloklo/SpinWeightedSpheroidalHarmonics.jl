# Spin-weighted spherical harmonics (direct sum, Jacobi recurrence, Chebyshev solve) and the
# Chebyshev solve of the spheroidal equation, all in the floating-point type of their arguments.
using LinearAlgebra
using ApproxFun
using SpecialFunctions
using AbstractTrees
include("binarytree.jl")

# Data of a node in the derivative tree built by _nth_derivative_spherical_harmonic_direct_eval:
# the term coeff * cos(theta/2)^ct2_power * sin(theta/2)^st2_power
struct ct2_st2
    coeff::Float64
    ct2_power::Int
    st2_power::Int
end

log_factorial(::Type{T}, n::Int) where {T} = n < 0 ? -T(Inf) : loggamma(T(n + 1))

function _log_summation_term_prefactor(::Type{T}, s::Int, l::Int, m::Int, r::Int) where {T}
    # Note that this does not include the (-1)^(l-r-s) factor
    # Check for negative arguments
    if (l-s-r) < 0 || (l-r+m) < 0 || (r+s-m) < 0
        # The whole thing is just 0
        return -T(Inf)
    end
    
    return begin 
        log_factorial(T, l-s) -
        log_factorial(T, l-s-r) - log_factorial(T, r) +
        log_factorial(T, l+s) -
        log_factorial(T, l-r+m) - log_factorial(T, r+s-m)
    end
end

function _summation_term_prefactors(::Type{T}, s::Int, l::Int, m::Int) where {T}
    log_prefactors = [_log_summation_term_prefactor(T, s, l, m, r) for r in 0:l-s]
    max_val, _ = findmax(log_prefactors)

    prefactor_signs = [(l-r-s) % 2 == 0 ? 1 : -1 for r in 0:l-s]
    log_prefactors = log_prefactors .- max_val # Now regularized
    prefactors = prefactor_signs .* exp.(log_prefactors)
    return prefactors, max_val
end

#=
sYlm(π/2, 0) and its θ-derivative there, in T. At the equator sin(θ/2) = cos(θ/2) = 1/√2, so
every term cos(θ/2)^a sin(θ/2)^b of the sum over r (a + b = 2l) is exactly 2^-l, and its
θ-derivative is (b - a)/2 times that, with (b - a)/2 = l - 2r - s + m: both sums are exact
integers over 2^l, summed exactly. The direct sum would lose digits to cancellation at large l.
The result is multiplied by the prefactor with at least 256 bits (and 64 guard bits beyond T),
whatever precision the caller has set.
=#
function _spherical_harmonic_at_pi_over_2(::Type{T}, s::Int, l::Int, m::Int; theta_derivative::Int=0) where {T}
    theta_derivative in (0, 1) || throw(ArgumentError("theta_derivative must be 0 or 1."))
    rmin = max(0, m-s)
    rmax = min(l-s, l+m)
    weight(r) = theta_derivative == 0 ? 1 : l - 2r - s + m
    terms = sum((iseven(l-r-s) ? 1 : -1) * weight(r) * binomial(big(l-s), r) * binomial(big(l+s), r+s-m)
        for r in rmin:rmax; init=big(0))
    prefactor = _swsh_prefactor(T, s, l, m)
    return setprecision(BigFloat, max(256, precision(T) + 64)) do
        T(BigFloat(prefactor) * (BigFloat(terms) / BigFloat(2)^l))
    end
end

#=
Chebyshev solve of the spheroidal ODE in x = cos(theta), split at the equator.

Near the poles S behaves like (1 - x)^a (1 + x)^b with a = |m + s|/2 and b = |m - s|/2,
which are half-integers when m ± s is odd, so we solve for the smooth F = S/((1 - x)^a (1 + x)^b)
instead. F satisfies a Jacobi-type equation,
    (1 - x^2) F'' + (2(b - a) - 2(a + b + 1) x) F' + R(x) F = 0,
and on each half we impose regularity at the pole (the equation there, where the F'' term drops out)
and one condition at the equator that fixes the scale. That condition weighs the value and the
derivative together, so that it does not vanish when the harmonic has a node at the equator.

`S_half` and `dS_half` are S and dS/dtheta at theta = π/2; `A` is the angular separation constant.
Return S(theta) as a Chebyshev series on [0, π].
=#
function _solve_spheroidal_harmonic_chebyshev(::Type{T}, s::Int, l::Int, m::Int, c, A, S_half, dS_half) where {T}
    a, b = T(abs(m + s)) / 2, T(abs(m - s)) / 2
    lmin = max(abs(m), abs(s))
    κ = l + 1 # typical size of d/dx, to weigh the value against the derivative
    # F and dF/dx at the equator, where (1 - x)^a (1 + x)^b = 1 and dS/dx = -dS/dtheta
    F0 = S_half
    F1 = -dS_half - (b - a) * S_half
    target = abs2(F0) + abs2(F1) / κ^2

    halves = map(((-one(T), zero(T)), (zero(T), one(T)))) do (lo, hi)
        dom = Chebyshev(lo..hi)
        x = Fun(dom)
        D = Derivative(dom)
        R = c^2 * x^2 - 2 * c * s * x + (A + s + s^2 - lmin * (lmin + 1))
        L = (1 - x^2) * D^2 + (2 * (b - a) - 2 * (a + b + 1) * x) * D + R
        E0 = Evaluation(dom, zero(T))
        normalization = conj(F0) * E0 + (conj(F1) / κ^2) * (E0 * D)
        pole = hi == 1 ? one(T) : -one(T)
        Ep = Evaluation(dom, pole)
        regularity = (2 * (b - a) - 2 * (a + b + 1) * pole) * (Ep * D) + R(pole) * Ep
        [normalization; regularity; L] \ [target, zero(target), zero(target)]
    end

    function S(θ)
        x = cos(θ)
        F = x < 0 ? halves[1](x) : halves[2](x)
        return (1 - x)^a * (1 + x)^b * F
    end
    return Fun(S, zero(T)..T(π))
end

function _solve_spherical_harmonic_chebyshev(::Type{T}, s::Int, l::Int, m::Int) where {T}
    # Evaluate sYlm(\pi/2, 0) using a numerically stable method
    Ypi2 = _spherical_harmonic_at_pi_over_2(T, s, l, m)
    dYpi2 = _spherical_harmonic_at_pi_over_2(T, s, l, m; theta_derivative=1)
    A = spin_weighted_spherical_eigenvalue(s, l, m)
    return _solve_spheroidal_harmonic_chebyshev(T, s, l, m, zero(T), A, Ypi2, dYpi2)
end

# Find the proper _theta in [0, π] and _phi to evaluate at, using the symmetry
# (theta, phi) -> (2π - theta, phi + π) of the harmonics.
function _fold_theta_phi(theta, phi)
    twopi = 2 * oftype(theta, π)
    _theta = mod(theta, twopi)
    _phi = phi
    _theta = _theta < 0 ? _theta + twopi : _theta
    if _theta > oftype(theta, π)
        _theta = twopi - _theta
        _phi += oftype(phi, π)
    end
    return _theta, _phi
end

function _nth_derivative_spheroidal_harmonic_chebyshev(chebyshev_S::Fun, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    _theta, _phi = _fold_theta_phi(theta, phi)

    # Note that the derivatives at the two boundary points might be inaccurate
    S_deriv = theta_derivative == 0 ? chebyshev_S(_theta) : differentiate(chebyshev_S, theta_derivative)(_theta)
    return S_deriv * cis(m * _phi) * (m * 1im)^phi_derivative
end

function _jacobi_polynomial_recurrence(n::Int, α::Int, β::Int, x::Real)
    if n < 0
        error("Jacobi polynomial order n must be non-negative")
    end

    xf = float(x)
    T = typeof(xf)
    αf = T(α)
    βf = T(β)

    n == 0 && return one(T)

    pnm1 = one(T)
    pn = 0.5 * ((2 + αf + βf) * xf + (αf - βf))
    n == 1 && return pn

    for k in 1:n-1
        kf = T(k)
        A = 2 * (kf + 1) * (kf + αf + βf + 1) * (2 * kf + αf + βf)
        B = (2 * kf + αf + βf + 1) * ((2 * kf + αf + βf + 2) * (2 * kf + αf + βf) * xf + αf^2 - βf^2)
        C = 2 * (kf + αf) * (kf + βf) * (2 * kf + αf + βf + 2)
        pnp1 = (B * pn - C * pnm1) / A
        pnm1, pn = pn, pnp1
    end
    return pn
end

#=
sin(θ/2)^a cos(θ/2)^b P_n^(α,β)(cos θ). A negative parameter α = -k makes P_n^(α,β) carry the
factor ((x - 1)/2)^k = (-sin²(θ/2))^k, which x = cos θ cannot resolve near θ = 0 (cos(1e-8) == 1):
the polynomial would be a cancellation down to roundoff. Szegő (4.22.2),
    C(n, k) P_n^(-k,β)(x) = C(n + β, k) ((x - 1)/2)^k P_{n-k}^(k,β)(x),
takes that factor out exactly in terms of sin(θ/2); likewise β = -k gives cos(θ/2)^(2k) through
P_n^(α,β)(-x) = (-1)^n P_n^(β,α)(x). (n ≥ k for every term of a harmonic.)
=#
function _jacobi_term(n::Int, α::Int, β::Int, a::Int, b::Int, st2, ct2, x)
    if α < 0
        k = -α
        ratio = prod(one(x) * (n + β - i) / (n - i) for i in 0:k-1)
        return (isodd(k) ? -ratio : ratio) * _jacobi_term(n - k, k, β, a + 2k, b, st2, ct2, x)
    elseif β < 0
        k = -β
        ratio = prod(one(x) * (n + α - i) / (n - i) for i in 0:k-1)
        return ratio * _jacobi_term(n - k, α, k, a, b + 2k, st2, ct2, x)
    end
    return st2^a * ct2^b * _jacobi_polynomial_recurrence(n, α, β, x)
end

function _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(s::Int, l::Int, m::Int, theta)
    m < 0 && error("Jacobi core expects m >= 0")

    α = m + s
    β = m - s
    n = l - m

    x = cos(theta)
    st2 = sin(theta / 2)
    ct2 = cos(theta / 2)

    return _swsh_prefactor(typeof(x), s, l, m) * _jacobi_term(n, α, β, α, β, st2, ct2, x)
end

function _differentiate_jacobi_terms(terms::Dict{NTuple{5, Int}, Float64})
    # Each term is represented as:
    # coeff * sin(theta/2)^a * cos(theta/2)^b * P_n^(p,q)(cos(theta))
    # with key (n, p, q, a, b).
    dterms = Dict{NTuple{5, Int}, Float64}()
    for (key, coeff) in terms
        n, p, q, a, b = key

        if a != 0
            key1 = (n, p, q, a - 1, b + 1)
            dterms[key1] = get(dterms, key1, 0.0) + coeff * (0.5 * a)
        end

        if b != 0
            key2 = (n, p, q, a + 1, b - 1)
            dterms[key2] = get(dterms, key2, 0.0) - coeff * (0.5 * b)
        end

        if n > 0
            key3 = (n - 1, p + 1, q + 1, a + 1, b + 1)
            dterms[key3] = get(dterms, key3, 0.0) - coeff * (n + p + q + 1)
        end
    end
    return dterms
end

# The theta-derivative of the m >= 0 Jacobi core as a list of terms coeff * sin(theta/2)^a * cos(theta/2)^b * P_n^(p,q)(cos theta),
# kept in the iteration order of the Dict that produced it (the evaluation sums in this order).
function _jacobi_theta_terms(s::Int, l::Int, m::Int, theta_derivative::Int)
    p = m + s
    q = m - s
    n = l - m
    terms = Dict{NTuple{5, Int}, Float64}((n, p, q, p, q) => 1.0)
    for _ in 1:theta_derivative
        terms = _differentiate_jacobi_terms(terms)
    end
    return collect(terms)
end

# The first and second theta-derivative term lists of a Jacobi-method harmonic, for the core parameters it evaluates
# (negative m is evaluated through the reflection (s, m) -> (-s, -m)).
function _jacobi_theta_derivative_terms(s::Int, l::Int, m::Int)
    ss, mm = m < 0 ? (-s, -m) : (s, m)
    return (_jacobi_theta_terms(ss, l, mm, 1), _jacobi_theta_terms(ss, l, mm, 2))
end

function _spin_weighted_spherical_harmonic_jacobi_theta_derivative_nonnegative_m(s::Int, l::Int, m::Int, theta, theta_derivative::Int;
        prepared_terms=nothing)
    m < 0 && error("Jacobi theta-derivative core expects m >= 0")
    theta_derivative < 0 && error("theta_derivative must be non-negative")

    terms = prepared_terms !== nothing && theta_derivative in (1, 2) ? prepared_terms[theta_derivative] :
        _jacobi_theta_terms(s, l, m, theta_derivative)

    x = cos(theta)
    st2 = sin(theta / 2)
    ct2 = cos(theta / 2)

    out = zero(x)
    for (key, coeff) in terms
        nterm, pterm, qterm, aterm, bterm = key
        out += coeff * _jacobi_term(nterm, pterm, qterm, aterm, bterm, st2, ct2, x)
    end

    return _swsh_prefactor(typeof(x), s, l, m) * out
end

function _spin_weighted_spherical_harmonic_jacobi_value(s::Int, l::Int, m::Int, theta, phi)
    # Fold theta into [0, π] using the symmetry (theta, phi) -> (2π - theta, phi + π) of the harmonics.
    _theta, _phi = _fold_theta_phi(theta, phi)
    T = typeof(_theta)

    # Exact values at the poles; a theta within 1e-14 of 0 or π (the tolerance scaled to the precision),
    # e.g. after the roundoff of the folding above, is treated as the pole.
    if isapprox(_theta, zero(T); atol=_tol(T, 1e-14))
        return (m == -s ? (-1)^s * sqrt((2*l+1)/(4*T(π))) : zero(T)) * cis(m * _phi)
    elseif isapprox(_theta, T(π); atol=_tol(T, 1e-14))
        return (m == s ? (-1)^l * sqrt((2*l+1)/(4*T(π))) : zero(T)) * cis(m * _phi)
    end

    if m >= 0
        return _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(s, l, m, _theta) * cis(m * _phi)
    end

    # Use the standard symmetry relation for negative m.
    reflected = _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(-s, l, -m, _theta) * cis((-m) * _phi)
    return (-1)^(s - m) * conj(reflected)
end

function _nth_derivative_spherical_harmonic_jacobi(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi;
        prepared_terms=nothing)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    # Fold theta into [0, π] using the symmetry (theta, phi) -> (2π - theta, phi + π) of the harmonics.
    _theta, _phi = _fold_theta_phi(theta, phi)

    if m < 0
        reflected = _nth_derivative_spherical_harmonic_jacobi(-s, l, -m, theta_derivative, phi_derivative, _theta, _phi;
            prepared_terms)
        return (-1)^(s - m) * conj(reflected)
    end

    if theta_derivative == 0
        y = _spin_weighted_spherical_harmonic_jacobi_value(s, l, m, _theta, _phi)
        return y * (m * 1im)^phi_derivative
    end

    # Near boundaries, use the direct expression to avoid catastrophic cancellation from negative powers.
    if isapprox(_theta, 0.0; atol=1e-3) || isapprox(_theta, π; atol=1e-3)
        return _nth_derivative_spherical_harmonic_direct_eval(s, l, m, theta_derivative, phi_derivative, _theta, _phi)
    end

    theta_part = _spin_weighted_spherical_harmonic_jacobi_theta_derivative_nonnegative_m(s, l, m, _theta, theta_derivative;
        prepared_terms)
    return theta_part * cis(m * _phi) * (m * 1im)^phi_derivative
end

#=
S and dS/dθ of a sum Σ_l w_l sYlm(θ, φ) over consecutive l of one (s, m), in one pass. The Jacobi
polynomials of every l come from one recurrence for the value and one for the derivative (their
parameters do not depend on l), run along with the sum, and the powers of sin(θ/2) and cos(θ/2) are
shared. The terms are those of _jacobi_term and _differentiate_jacobi_terms, using
dP_n^(α,β)(cos θ)/dθ = -(n + α + β + 1)/2 sin θ P_{n-1}^(α+1,β+1)(cos θ): a negative parameter is taken
out by the same identity, so no power is negative and the value and the derivative stay accurate near
the poles. Negative m goes through the reflection (s, m) -> (-s, -m), as in
_spin_weighted_spherical_harmonic_jacobi_value; the core is then real, so the reflection only changes
the overall sign and the phase. `weight(l)` is w_l, or `nothing` for a harmonic left out of the sum;
theta must lie in (0, π), away from the poles (see _near_pole). Nothing is allocated in Float64.
=#
function _jacobi_value_and_theta_derivative(weight, s::Int, m::Int, lmin::Int, lmax::Int, theta, phi)
    sc, mc = m < 0 ? (-s, -m) : (s, m)
    α, β = mc + sc, mc - sc
    x = cos(theta)
    st2 = sin(theta / 2)
    ct2 = cos(theta / 2)
    T = typeof(x)
    α0, β0, da0, db0, k0 = _jacobi_shift(α, β)
    α1, β1, da1, db1, k1 = _jacobi_shift(α + 1, β + 1)
    # the powers of sin(θ/2) and cos(θ/2) of the value and of the three terms of the derivative
    sv, cv = st2^(α + da0), ct2^(β + db0)
    sa, ca = α != 0 ? st2^(α - 1 + da0) : zero(T), ct2^(β + 1 + db0)
    sb, cb = st2^(α + 1 + da0), β != 0 ? ct2^(β - 1 + db0) : zero(T)
    sd, cd = st2^(α + 1 + da1), ct2^(β + 1 + db1)
    # P_j0 of the value and P_j1 of the derivative, with their predecessors
    j0, p0, q0 = 0, one(T), zero(T)
    j1, p1, q1 = 0, one(T), zero(T)
    value = derivative = zero(T) * zero(weight(lmin) === nothing ? zero(T) : weight(lmin))
    for l in lmin:lmax
        w = weight(l)
        w === nothing && continue
        n = l - mc
        while j0 < n - k0
            j0, p0, q0 = j0 + 1, _jacobi_next(j0, p0, q0, α0, β0, x), p0
        end
        Q = _jacobi_shift_ratio(n, α, β, x) * p0
        d = zero(T)
        α != 0 && (d += (T(α) / 2) * sa * ca * Q)
        β != 0 && (d -= (T(β) / 2) * sb * cb * Q)
        if n > 0
            while j1 < n - 1 - k1
                j1, p1, q1 = j1 + 1, _jacobi_next(j1, p1, q1, α1, β1, x), p1
            end
            d -= (n + α + β + 1) * sd * cd * _jacobi_shift_ratio(n - 1, α + 1, β + 1, x) * p1
        end
        prefactor = _swsh_prefactor(T, sc, l, mc)
        value += w * (prefactor * (sv * cv * Q))
        derivative += w * (prefactor * d)
    end
    phase = (m < 0 ? (-1)^(s - m) : 1) * cis(m * phi)
    return value * phase, derivative * phase
end

# P_{j+1}^(α,β)(x) from P_j = pj and P_{j-1} = pjm1, by the recurrence of _jacobi_polynomial_recurrence
@inline function _jacobi_next(j::Int, pj, pjm1, α::Int, β::Int, x)
    T = typeof(x)
    αf, βf = T(α), T(β)
    j == 0 && return 0.5 * ((2 + αf + βf) * x + (αf - βf))
    kf = T(j)
    A = 2 * (kf + 1) * (kf + αf + βf + 1) * (2 * kf + αf + βf)
    B = (2 * kf + αf + βf + 1) * ((2 * kf + αf + βf + 2) * (2 * kf + αf + βf) * x + αf^2 - βf^2)
    C = 2 * (kf + αf) * (kf + βf) * (2 * kf + αf + βf + 2)
    return (B * pj - C * pjm1) / A
end

# The identity of _jacobi_term: P_n^(α,β) = ratio(n) sin(θ/2)^da cos(θ/2)^db P_{n-k}^(α',β'), with
# (α', β', da, db, k)
_jacobi_shift(α::Int, β::Int) = α < 0 ? (-α, β, -2α, 0, -α) : β < 0 ? (α, -β, 0, -2β, -β) : (α, β, 0, 0, 0)

function _jacobi_shift_ratio(n::Int, α::Int, β::Int, x)
    ratio = one(x)
    if α < 0
        for i in 0:-α-1
            ratio *= one(x) * (n + β - i) / (n - i)
        end
        return isodd(α) ? -ratio : ratio
    elseif β < 0
        for i in 0:-β-1
            ratio *= one(x) * (n + α - i) / (n - i)
        end
    end
    return ratio
end

# Within the pole tolerance of _spin_weighted_spherical_harmonic_jacobi_value, which takes the exact
# pole value there
_near_pole(theta) = (T = typeof(theta); isapprox(theta, zero(T); atol=_tol(T, 1e-14)) ||
    isapprox(theta, T(π); atol=_tol(T, 1e-14)))

# The terms coeff * cos(θ/2)^a * sin(θ/2)^b of the theta_derivative-th derivative of
# cos(θ/2)^a0 * sin(θ/2)^b0, in the order the derivative tree lists its leaves, zero terms dropped
# (the coefficients are exact: products of ±1/2 and integers)
function _direct_derivative_terms(a0::Int, b0::Int, theta_derivative::Int)
    root = BinaryNode(ct2_st2(1, a0, b0)) # root of the tree
    # building the binary tree
    for j in 1:theta_derivative
        #=
            Each derivative wrt theta will add two terms
            one with
            1/2 beta ct2^{alpha+1} st2^{beta-1}
            and one with
            -1/2 alpha ct2^{alpha-1} st2^{beta+1}
        =#
        # now find the appropriate parent
        # traverse the current tree looking for leaves to add new nodes
        for leaf in Leaves(root)
            leftchild(ct2_st2(leaf.data.coeff * 0.5 * leaf.data.st2_power, leaf.data.ct2_power+1, leaf.data.st2_power-1), leaf)
            rightchild(ct2_st2(leaf.data.coeff * -0.5 * leaf.data.ct2_power, leaf.data.ct2_power-1, leaf.data.st2_power+1), leaf)
        end
    end
    return [(leaf.data.coeff, leaf.data.ct2_power, leaf.data.st2_power) for leaf in Leaves(root) if leaf.data.coeff != 0]
end

#=
What the direct sum needs besides the angle, prepared once with the harmonic in its type T: the
normalized prefactors of the sum over r, their log normalization, the overall prefactor, and the
terms of the value and of its first two theta-derivatives for every r.
=#
struct DirectSumTerms{T}
    prefactors::Vector{T}
    log_normalization::T
    swsh_prefactor::T
    terms::NTuple{3, Vector{Vector{Tuple{Float64, Int, Int}}}}
end

function DirectSumTerms(::Type{T}, s::Int, l::Int, m::Int) where {T}
    prefactors, log_normalization = _summation_term_prefactors(T, s, l, m)
    rs = max(0, m-s):min(l-s, l+m)
    terms = ntuple(d -> [_direct_derivative_terms(2*r+s-m, 2*l-2*r-s+m, d - 1) for r in rs], 3)
    return DirectSumTerms{T}(prefactors, log_normalization, _swsh_prefactor(T, s, l, m), terms)
end

function _nth_derivative_spherical_harmonic_direct_eval(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi;
        prepared=nothing)
    # Outside [0, π] the harmonic is evaluated at the mirrored point, as in the other methods
    theta, phi = _fold_theta_phi(theta, phi)
    ct2 = cos(theta/2)
    st2 = sin(theta/2)

    T = typeof(ct2)
    # Prepared terms are used when they are in the type of the angle and hold this derivative
    usable = prepared isa DirectSumTerms{T} && theta_derivative <= 2
    summation_term_prefactors, log_normalization_const = usable ?
        (prepared.prefactors, prepared.log_normalization) : _summation_term_prefactors(T, s, l, m)
    _sum = zero(T)
    for (k, r) in enumerate(max(0, m-s):min(l-s, l+m))
        _rsum = zero(T)
        terms = usable ? prepared.terms[theta_derivative + 1][k] :
            _direct_derivative_terms(2*r+s-m, 2*l-2*r-s+m, theta_derivative)
        for (coeff, ct2_power, st2_power) in terms
            _rsum += coeff * ct2^ct2_power * st2^st2_power
        end
        _rsum *= summation_term_prefactors[r+1]
        _sum += _rsum
    end

    # Check if _sum is zero before taking log
    if _sum == 0
        return zero(T)
    end

    swsh_pref = usable ? prepared.swsh_prefactor : _swsh_prefactor(T, s, l, m)
    log_amp = log_normalization_const + log(abs(_sum)) + log(abs(swsh_pref))
    sign_total = sign(_sum) * sign(swsh_pref)
    amp = exp(log_amp)
    return amp * sign_total * cis(m*phi) * (m*1im)^phi_derivative

end

function _swsh_prefactor(::Type{T}, s::Int, l::Int, m::Int) where {T}
    #=
        This is consistent with the expression on Wikipedia
        and in the Black Hole Perturbation Toolkit.
    =#

    common_factor = (-1)^m * sqrt((2*l+1)/(4*T(π)))
    if abs(s) == abs(m)
        return common_factor
    elseif s > m
        delta = s - m # which is a positive integer
        out = one(T)
        for i in 0:1:delta-1
            j = delta - i
            out *= (T(l-m-i)/T(l+m+j))
        end
        return common_factor * sqrt(out)
    else
        # in this case s < m
        delta = m - s # which is a positive integer
        out = one(T)
        for i in 0:1:delta-1
            j = delta - i
            out *= (T(l+s+j)/T(l-s-i))
        end
        return common_factor * sqrt(out)
    end
end
