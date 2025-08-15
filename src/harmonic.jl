using LinearAlgebra
using ApproxFun
using TaylorSeries
using LogarithmicNumbers
using AbstractTrees
include("binarytree.jl")

# A tree node
struct ct2_st2
    coeff::Float64
    ct2_power::Int
    st2_power::Int
end

function log_factorial(n::Int)
    if n == 0
        return 0
    elseif n < 0
        return -Inf
    else
        return sum(log(k) for k in 1:n)
    end
end

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
    two_to_the_pow_2 = ULogarithmic(BigInt(2))^(ULogarithmic(BigInt(l)))
    _swsh_prefactor(s, l, m)*float(sum([prefactor_sign(s,l,r)*_summation_term_prefactor(s,l,m,r)/two_to_the_pow_2 for r in rmin:1:rmax]))
end

function _solve_spherical_harmonic_chebyshev(s::Int, l::Int, m::Int)
    # Evaluate sYlm(0,0) using the direct method for consistency
    Y0 = real(_nth_derivative_spherical_harmonic_direct_eval(s, l, m, 0, 0, 0, 0))
    # Evaluate sYlm(\pi/2, 0) using a numerically stable method
    Ypi2 = Float64(spin_weighted_spherical_harmonic_at_pi_over_2(s, l, m))
    Ypi = 0.0 # Always 0

    # Split the domain into two parts -- one from theta = π to π/2 and from π/2 to 0
    # NOTE x=cos(theta), x = -1 when \theta is \pi and x = 1 when \theta is 0

    # Solve in the first domain
    a, b = -1, 0
    dom1 = a..b

    # Define the differential operator
    x = Fun(dom1)
    D = Derivative(dom1)
    L = (1 - x^2)*D^2 - 2*x*D + (s + l*(l+1) - s*(s+1) - (m + s*x)^2/(1-x^2))

    bvals = [Ypi, Ypi2] # Boundary values
    u = [Dirichlet(dom1); L] \ [bvals, 0]

    # Solve in the second domain
    dom2 = 0..1
    x = Fun(dom2)
    D = Derivative(dom2)
    L = (1 - x^2)*D^2 - 2*x*D + (s + l*(l+1) - s*(s+1) - (m + s*x)^2/(1-x^2))

    bvals = [Ypi2, Y0] # Boundary values
    v = [Dirichlet(dom2); L] \ [bvals, 0]

    Y(θ) = θ > π/2 ? u(cos(θ)) : v(cos(θ))
    return Y
end

function _nth_derivative_spherical_harmonic_chebyshev(chebyshev_Y::Function, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    factorial(theta_derivative)*getcoeff(taylor_expand(chebyshev_Y, theta, order=theta_derivative), theta_derivative) * cis(m*phi) * (m*1im)^phi_derivative
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
    exp(log_normalization_const) * _swsh_prefactor(s, l, m) * _sum * cis(m*phi) * (m*1im)^phi_derivative
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
