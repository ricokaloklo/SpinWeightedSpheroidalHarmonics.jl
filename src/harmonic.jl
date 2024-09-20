using LinearAlgebra
using ApproxFun
using TaylorSeries
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

function _nth_derivative_spherical_harmonic(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi; method="auto")
    if method == "auto"
        _method = l >= 30 ? "chebyshev" : "direct"
        return _nth_derivative_spherical_harmonic(s, l, m, theta_derivative, phi_derivative, theta, phi; method=_method)
    elseif method == "direct"
        return _nth_derivative_spherical_harmonic_direct_eval(s, l, m, theta_derivative, phi_derivative, theta, phi)
    elseif method == "chebyshev"
        return _nth_derivative_spherical_harmonic_chebyshev(s, l, m, theta_derivative, phi_derivative, theta, phi)
    else
        error("Does not understand method $method")
    end
end

function _nth_derivative_spherical_harmonic_chebyshev(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    # Evaluate sYlm(0,0) using the direct method for consistency
    Y0 = real(_nth_derivative_spherical_harmonic_direct_eval(s, l, m, 0, 0, 0, 0))

    # Define the domain
    # NOTE x=cos(theta), x = -1 when \theta is \pi and x = 1 when \theta is 0
    a, b = -1, 1;
    dom = a..b;

    # Define the differential operator
    x = Fun(dom);
    D = Derivative(dom);
    L = (1 - x^2)*D^2 - 2*x*D + (s + l*(l+1) - s*(s+1) - (m + s*x)^2/(1-x^2));

    bvals = [0, Y0] # Boundary values
    u = [Dirichlet(dom); L] \ [bvals, 0];
    Y(θ) = u(cos(θ))

    factorial(theta_derivative)*getcoeff(taylor_expand(Y, theta, order=theta_derivative), theta_derivative) * cis(m*phi) * (m*1im)^phi_derivative
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
            _rsum += leaf.data.coeff * ct2^(leaf.data.ct2_power) * st2^(leaf.data.st2_power)
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
