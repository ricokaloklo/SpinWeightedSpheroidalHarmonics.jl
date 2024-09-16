function _determine_matrix_size_N(s::Int, l::Int, m::Int)
    #=
    Determine a suitable value of N for the spectral decomposition

    The value of N calculated here is essentially lmax for
    the spectral decomposition. Then we apply a 'buffer' of 10
    =#
    N = l - max(abs(m), abs(s)) + 1
    return N + 10
end
