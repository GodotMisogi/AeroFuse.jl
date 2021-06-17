# Old (Grandpa) approach
function embedding((m, n), block :: Matrix{<: Real}, c_index)
    mat = zeros(m, n)
    mat[c_index .+ CartesianIndices(block) .- CartesianIndex(1,1)] .= block

    mat
end

beam_stiffness_matrix(n, E, Iz, L) = sum(embedding((2(n+1), 2(n+1)), Iz_coeffs(E, Iz, L), c_inds) for c_inds in CartesianIndex.(1:2:2n, 1:2:2n))


    K = zeros(6(n+1), 6(n+1))
    J_coeff_size = CartesianIndices((2, 2))
    I_coeff_size = CartesianIndices((4, 4))

    Iy_indices = CartesianIndex(0, 0)
    Iz_indices = CartesianIndex(2(n+1), 2(n+1))
    A_indices  = CartesianIndex(4(n+1), 4(n+1))
    J_indices  = CartesianIndex(5(n+1), 5(n+1))

    for (c_index, E, Iz, Iy, L) ∈ zip(CartesianIndex.(0:2:2(n-1), 0:2:2(n-1)), Es, Izs, Iys, Ls)
        K[c_index + Iy_indices .+ I_coeff_size] .+= Iy_coeffs(E, Iy, L)
        K[c_index + Iz_indices .+ I_coeff_size] .+= Iz_coeffs(E, Iz, L)
    end

    for (c_index, E, A, G, J, L) ∈ zip(CartesianIndex.(0:1:n-1, 0:1:n-1), Es, As, Gs, Js, Ls)
        K[c_index + A_indices .+ J_coeff_size] .+= J_coeffs(E, A, L)
        K[c_index + J_indices .+ J_coeff_size] .+= J_coeffs(G, J, L)
    end

    sparse(K)
end 