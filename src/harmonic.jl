using LinearAlgebra
using ApproxFun
using SpecialFunctions
using LogarithmicNumbers
using AbstractTrees
include("binarytree.jl")

# A tree node
struct ct2_st2
    coeff::Float64
    ct2_power::Int
    st2_power::Int
end

log_factorial(n::Int) = n < 0 ? -Inf : loggamma(n + 1)

function _log_summation_term_prefactor(s::Int, l::Int, m::Int, r::Int)
    # Note that this does not include the (-1)^(l-r-s) factor
    # Check for negative arguments
    if (l-s-r) < 0 || (l-r+m) < 0 || (r+s-m) < 0
        # The whole thing is just 0
        return -Inf
    end
    
    return begin 
        log_factorial(l-s) -
        log_factorial(l-s-r) - log_factorial(r) +
        log_factorial(l+s) -
        log_factorial(l-r+m) - log_factorial(r+s-m)
    end
end

function _summation_term_prefactors(s::Int, l::Int, m::Int)
    log_prefactors = [_log_summation_term_prefactor(s, l, m, r) for r in 0:l-s]
    max_val, _ = findmax(log_prefactors)

    prefactor_signs = [(l-r-s) % 2 == 0 ? 1 : -1 for r in 0:l-s]
    log_prefactors = log_prefactors .- max_val # Now regularized
    prefactors = prefactor_signs .* exp.(log_prefactors)
    return prefactors, max_val
end

function spin_weighted_spherical_harmonic_at_pi_over_2(s::Int, l::Int, m::Int)
    # NOTE These functions might look familiar to you
    # Indeed they are the same as above, but using LogarithmicNumbers with arbitrary precisions
    # NOTE They are not fast
    function bigfactorial(n)
        if n == 0
            return ULogarithmic(BigInt(1))
        end

        return prod(ULogarithmic(BigInt(k)) for k in 1:n)
    end

    function _summation_term_prefactor(s::Int, l::Int, m::Int, r::Int)
        # Note that this does not include the (-1)^(l-r-s) factor
        # Check for negative arguments
        if (l-s-r) < 0 || (l-r+m) < 0 || (r+s-m) < 0
            # The whole thing is just 0
            return 0
        end

        return (bigfactorial(l-s)/(bigfactorial(l-s-r)*bigfactorial(r))) * (bigfactorial(l+s)/(bigfactorial(l-r+m)*bigfactorial(r+s-m)))
    end

    rmin = max(0, m-s)
    rmax = min(l-s, l+m)

    prefactor_sign(s, l, r) = (l-r-s) % 2 == 0 ? 1 : -1
    two_to_the_pow_l = ULogarithmic(BigInt(2))^(ULogarithmic(BigInt(l)))
    _swsh_prefactor(s, l, m)*float(sum([prefactor_sign(s,l,r)*_summation_term_prefactor(s,l,m,r)/two_to_the_pow_l for r in rmin:1:rmax]))
end

function _solve_spheroidal_harmonic_chebyshev(s::Int, m::Int, c, λ, S0, Spi2, Spi)
    # Split the domain into two parts -- one from theta = π to π/2 and from π/2 to 0
    # NOTE x=cos(theta), x = -1 when \theta is \pi and x = 1 when \theta is 0

    # Solve in the first domain
    dom1 = -1..0

    # Define the differential operator
    x = Fun(dom1)
    D = Derivative(dom1)
    L = (1 - x^2) * (1 - x^2) * D^2 - 2 * x * (1 - x^2) * D + (((c * x)^2 - 2 * c * s * x + s + λ - c^2 + 2 * m * c) * (1 - x^2) - (m + s * x)^2)

    rhs_zero = zero(S0 + Spi2 + Spi)
    bvals = [Spi, Spi2] # Boundary values
    u = [Dirichlet(dom1); L] \ [bvals, rhs_zero]

    # Solve in the second domain
    dom2 = 0..1
    x = Fun(dom2)
    D = Derivative(dom2)
    L = (1 - x^2) * (1 - x^2) * D^2 - 2 * x * (1 - x^2) * D + (((c * x)^2 - 2 * c * s * x + s + λ - c^2 + 2 * m * c) * (1 - x^2) - (m + s * x)^2)

    bvals = [Spi2, S0] # Boundary values
    v = [Dirichlet(dom2); L] \ [bvals, rhs_zero]

    S(θ) = θ > π / 2 ? u(cos(θ)) : v(cos(θ))
    return S
end

function _solve_spherical_harmonic_chebyshev(s::Int, l::Int, m::Int)
    # Evaluate sYlm(0,0) exactly
    Y0 = m == -s ? (-1)^s * sqrt((2*l+1)/(4π)) : 0.0
    # Evaluate sYlm(\pi/2, 0) using a numerically stable method
    Ypi2 = Float64(spin_weighted_spherical_harmonic_at_pi_over_2(s, l, m))
    # Evaluate sYlm(\pi,0) exactly
    Ypi = m == s ? (-1)^l * sqrt((2*l+1)/(4π)) : 0.0

    return _solve_spheroidal_harmonic_chebyshev(s, m, 0.0, spin_weighted_spherical_eigenvalue(s, l, m), Y0, Ypi2, Ypi)
end

function _nth_derivative_spheroidal_harmonic_chebyshev(chebyshev_S::Fun, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    # Find the proper _theta in [0, π] and _phi in [0, 2π) to evaluate.
    _theta = mod(theta, 2π)
    _phi = phi
    _theta = _theta < 0 ? _theta + 2π : _theta
    if _theta > π
        _theta = 2π - _theta
        _phi += π
    end

    # Note that the derivatives at the two boundary points might be inaccurate
    S_deriv = theta_derivative == 0 ? chebyshev_S(_theta) : differentiate(chebyshev_S, theta_derivative)(_theta)
    return S_deriv * cis(m * _phi) * (m * 1im)^phi_derivative
end

function _nth_derivative_spheroidal_harmonic_chebyshev(chebyshev_S::Function, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    S_fun = Fun(chebyshev_S, 0..π)
    return _nth_derivative_spheroidal_harmonic_chebyshev(S_fun, m, theta_derivative, phi_derivative, theta, phi)
end

function _jacobi_polynomial_recurrence(n::Int, α::Int, β::Int, x::Real)
    if n < 0
        error("Jacobi polynomial order n must be non-negative")
    end

    αf = float(α)
    βf = float(β)
    xf = float(x)

    n == 0 && return 1.0

    pnm1 = 1.0
    pn = 0.5 * ((2 + αf + βf) * xf + (αf - βf))
    n == 1 && return pn

    for k in 1:n-1
        kf = float(k)
        A = 2 * (kf + 1) * (kf + αf + βf + 1) * (2 * kf + αf + βf)
        B = (2 * kf + αf + βf + 1) * ((2 * kf + αf + βf + 2) * (2 * kf + αf + βf) * xf + αf^2 - βf^2)
        C = 2 * (kf + αf) * (kf + βf) * (2 * kf + αf + βf + 2)
        pnp1 = (B * pn - C * pnm1) / A
        pnm1, pn = pn, pnp1
    end
    return pn
end

function _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(s::Int, l::Int, m::Int, theta)
    m < 0 && error("Jacobi core expects m >= 0")

    α = m + s
    β = m - s
    n = l - m

    x = cos(theta)
    st2 = sin(theta / 2)
    ct2 = cos(theta / 2)
    P = _jacobi_polynomial_recurrence(n, α, β, x)

    return _swsh_prefactor(s, l, m) * st2^α * ct2^β * P
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

function _spin_weighted_spherical_harmonic_jacobi_theta_derivative_nonnegative_m(s::Int, l::Int, m::Int, theta, theta_derivative::Int)
    m < 0 && error("Jacobi theta-derivative core expects m >= 0")
    theta_derivative < 0 && error("theta_derivative must be non-negative")

    p = m + s
    q = m - s
    n = l - m

    terms = Dict{NTuple{5, Int}, Float64}((n, p, q, p, q) => 1.0)
    for _ in 1:theta_derivative
        terms = _differentiate_jacobi_terms(terms)
    end

    x = cos(theta)
    st2 = sin(theta / 2)
    ct2 = cos(theta / 2)

    out = 0.0
    for (key, coeff) in terms
        nterm, pterm, qterm, aterm, bterm = key
        out += coeff * st2^aterm * ct2^bterm * _jacobi_polynomial_recurrence(nterm, pterm, qterm, x)
    end

    return _swsh_prefactor(s, l, m) * out
end

function _spin_weighted_spherical_harmonic_jacobi_value(s::Int, l::Int, m::Int, theta, phi)
    # Map angles to the principal interval to avoid branch issues from integer powers.
    _theta = mod(theta, 2π)
    _phi = phi
    _theta = _theta < 0 ? _theta + 2π : _theta
    if _theta > π
        _theta = 2π - _theta
        _phi += π
    end

    # Exact boundary values.
    if isapprox(_theta, 0.0; atol=1e-14)
        return (m == -s ? (-1)^s * sqrt((2*l+1)/(4π)) : 0.0) * cis(m * _phi)
    elseif isapprox(_theta, π; atol=1e-14)
        return (m == s ? (-1)^l * sqrt((2*l+1)/(4π)) : 0.0) * cis(m * _phi)
    end

    if m >= 0
        return _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(s, l, m, _theta) * cis(m * _phi)
    end

    # Use the standard symmetry relation for negative m.
    reflected = _spin_weighted_spherical_harmonic_jacobi_core_nonnegative_m(-s, l, -m, _theta) * cis((-m) * _phi)
    return (-1)^(s - m) * conj(reflected)
end

function _nth_derivative_spherical_harmonic_jacobi(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    theta_derivative < 0 && error("theta_derivative must be non-negative")
    phi_derivative < 0 && error("phi_derivative must be non-negative")

    # Find the proper _theta in [0, π] and _phi in [0, 2π) to evaluate.
    _theta = mod(theta, 2π)
    _phi = phi
    _theta = _theta < 0 ? _theta + 2π : _theta
    if _theta > π
        _theta = 2π - _theta
        _phi += π
    end

    if m < 0
        reflected = _nth_derivative_spherical_harmonic_jacobi(-s, l, -m, theta_derivative, phi_derivative, _theta, _phi)
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

    theta_part = _spin_weighted_spherical_harmonic_jacobi_theta_derivative_nonnegative_m(s, l, m, _theta, theta_derivative)
    return theta_part * cis(m * _phi) * (m * 1im)^phi_derivative
end

function _nth_derivative_spherical_harmonic_direct_eval(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    ct2 = cos(theta/2)
    st2 = sin(theta/2)

    _sum = 0.0
    summation_term_prefactors, log_normalization_const = _summation_term_prefactors(s, l, m)
    for r in max(0, m-s):min(l-s, l+m)
        _rsum = 0.0
        root = BinaryNode(ct2_st2(1, 2*r+s-m, 2*l-2*r-s+m)) # root of the tree
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
        # now traverse the final tree
        for leaf in Leaves(root)
            if leaf.data.coeff != 0
                _rsum += leaf.data.coeff * ct2^(leaf.data.ct2_power) * st2^(leaf.data.st2_power)
            end
        end
        _rsum *= summation_term_prefactors[r+1]
        _sum += _rsum
    end

    # Check if _sum is zero before taking log
    if _sum == 0.0
        return 0.0
    end

    swsh_pref = _swsh_prefactor(s, l, m)
    log_amp = log_normalization_const + log(abs(_sum)) + log(abs(swsh_pref))
    sign_total = sign(_sum) * sign(swsh_pref)
    amp = exp(log_amp)
    return amp * sign_total * cis(m*phi) * (m*1im)^phi_derivative

end

function _swsh_prefactor(s::Int, l::Int, m::Int)
    #=
        This is consistent with the expression in wikipedia,
        as well as BHPerturbationToolkit
    =#

    # Implement explicit expression here
    common_factor = (-1)^m * sqrt((2*l+1)/(4*pi))
    if abs(s) == abs(m)
        return common_factor
    elseif s > m
        delta = s - m # which is a positive integer
        out = 1
        for i in 0:1:delta-1
            j = delta - i
            out *= ((l-m-i)/(l+m+j))
        end
        return common_factor * sqrt(out)
    else
        # in this case s < m
        delta = m - s # which is a positive integer
        out = 1
        for i in 0:1:delta-1
            j = delta - i
            out *= ((l+s+j)/(l-s-i))
        end
        return common_factor * sqrt(out)
    end
end
