# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(r, horseshoe, normal, V_hat, symmetry)

Computes the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
function influence_coefficient(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, normal :: SVector{3, <: Real}, V_hat :: SVector{3, <: Real}, symmetry :: Bool)
    if symmetry
        @timeit "Velocity" col_vel = velocity(r, horseshoe, 1., V_hat)
        @timeit "Ref Velocity" ref_vel = (reflect_xz ∘ velocity)(reflect_xz(r), horseshoe, 1., V_hat)

        @timeit "Dotting" dot(col_vel + ref_vel, normal)
    else
        @timeit "Dotting" dot(velocity(r, horseshoe, 1., V_hat), normal)
    end
end

"""
    influence_matrix(colpoints, normals, horseshoes, V_hat, symmetry)

Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, and a unit vector representing the freestream.
"""
influence_matrix(colpoints, normals, horseshoes :: AbstractVector{<: Horseshoe}, V_hat :: SVector{3, <: Real}, symmetry :: Bool) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(r_i, horsie_j, n_i, V_hat, symmetry) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]

"""
    boundary_condition(velocities, normals)

Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)

# Pre-allocation attempts
#==========================================================================================#

function influence_matrix!(AIC :: AbstractMatrix{<: Real}, colpoints :: AbstractArray{<: SVector{3, <: Real}}, normals :: AbstractArray{<: SVector{3, <: Real}}, horseshoes :: AbstractVector{<: Horseshoe}, V_hat :: SVector{3, <: Real})
    for i ∈ 1:length(colpoints)
        for j ∈ 1:length(horseshoes)
            @timeit "a" a = r1(colpoints[i], horseshoes[j])
            @timeit "b" b = r2(colpoints[i], horseshoes[j])
            @timeit "Allocating?" AIC[i,j] = dot(total_horseshoe_velocity(a, b, 1., V_hat), normals[i])
        end
    end
    AIC
end

function boundary_condition!(boco, velocities, normals)
    @. boco = dot(velocities, normals)
    boco
end