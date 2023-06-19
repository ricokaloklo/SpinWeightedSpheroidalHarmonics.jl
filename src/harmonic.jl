using AbstractTrees
include("binarytree.jl")

# A tree node
struct ct2_st2
    coeff::Float64
    ct2_power::Int
    st2_power::Int
end

function _summation_term_prefactor(s::Int, l::Int, m::Int, r::Int)
    binomial(l-s, r)*binomial(l+s, r+s-m)*(-1)^(l-r-s)
end

function _nth_derivative_spherical_harmonic(s::Int, l::Int, m::Int, theta_derivative::Int, phi_derivative::Int, theta, phi)
    ct2 = cos(theta/2)
    st2 = sin(theta/2)

    _sum = 0.0
    for r in 0:l-s
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
        _rsum *= _summation_term_prefactor(s, l, m, r)
        _sum += _rsum
    end
    _swsh_prefactor(s, l, m) * _sum * cis(m*phi) * (m*1im)^phi_derivative
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
        diff = s - m # which is a positive integer
        out = 1
        for i in 0:1:diff-1
            out *= (l-m-i)
        end
        for j in diff:-1:1
            out /= (l+m+j)
        end
        return common_factor * sqrt(out)
    else
        # in this case s < m
        diff = m - s # which is a positive integer
        out = 1
        for i in diff:-1:1
            out *= (l+s+i)
        end
        for j in 0:1:diff-1
            out /= (l-s-j)
        end
        return common_factor * sqrt(out)
    end
    # should throw an error here because the function should have exited before this
end

@doc raw"""
    spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)

Compute the spin-weighted spherical harmonic with spin weight `s`, harmonic index `l`, azimuthal index `m`, and coordinates `theta` and `phi`.

The optional arguments `theta_derivative` and `phi_derivative` specify the order of partial derivatives to take with respect to `theta` and `phi`, respectively.
"""
function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)
    _nth_derivative_spherical_harmonic(s, l, m, theta_derivative, phi_derivative, theta, phi)
end