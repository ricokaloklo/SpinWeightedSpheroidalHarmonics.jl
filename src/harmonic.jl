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
    (-1)^m * sqrt(Complex((factorial(big(l+m))*factorial(big(l-m))*(2*l+1))/(4*pi*factorial(big(l+s))*factorial(big(l-s)))))
end

function spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int, theta, phi)
    _nth_derivative_spherical_harmonic(s, l, m, 0, 0, theta, phi)
end