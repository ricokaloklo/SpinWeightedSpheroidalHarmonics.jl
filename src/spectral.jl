using LinearAlgebra

function Fslm(s::Int, l::Int, m::Int)
    sqrt( ( (l+1)^2 - m^2 )/( (2*l+3)*(2*l+1) ) ) * sqrt( ( (l+1)^2 - s^2 )/( (l+1)^2 ) )
end

function Gslm(s::Int, l::Int, m::Int)
    l == 0 ? 0 : sqrt( ( l^2 - m^2 )/( 4*l^2 - 1 ) ) * sqrt( ( l^2 - s^2 )/( l^2 ) )
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

function angular_sep_const(c, s::Int, l::Int, m::Int, N::Int=10)
    all_l_in_matrix = construct_all_l_in_matrix(s, m, N)
    spectral_matrix = construct_spectral_matrix(c, s, m, N)

    # Solve the eigenvalue problem
    eigenvalues = eigvals(spectral_matrix)
    # The returned eigenvalues are already sorted in ascending order
    eigenvalues[indexin(l, all_l_in_matrix)[1]]
end

function spectral_coefficients(c, s::Int, l::Int, m::Int, N::Int=10)
    all_l_in_matrix = construct_all_l_in_matrix(s, m, N)
    spectral_matrix = construct_spectral_matrix(c, s, m, N)

    # Solve the eigenvalue problem
    eigenvectors = eigvecs(spectral_matrix)
    # The returned eigenvectors are already normalized
    v = eigenvectors[:,indexin(l, all_l_in_matrix)[1]]
    # Note that the eigenvectors are determined only up to a multiplicative factor
    # Make sure that when l'=l, it is a positive real number
    v /= v[indexin(l, all_l_in_matrix)[1]]
    # Re-normalize the vector
    return v/sqrt(dot(v,v))
end

function Teukolsky_lambda_const(c, s::Int, l::Int, m::Int, N::Int=10)
    angular_sep_const(c, s, l, m, N) + c^2 - 2*m*c
end