# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(r, horseshoe, panel_normal, V_hat, symmetry)

Computes the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
function influence_coefficient(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, panel_normal :: SVector{3, <: Real}, V_hat :: SVector{3, <: Real}, symmetry :: Bool)
    if symmetry
        col_vel = velocity(r, horseshoe, 1., V_hat)
        mir_vel = (reflect_xz ∘ velocity)(reflect_xz(r), horseshoe, 1., V_hat)

        dot(col_vel + mir_vel, panel_normal)
    else
        dot(velocity(r, horseshoe, 1., V_hat), panel_normal)
    end
end


"""
    influence_matrix(colpoints, normals, horseshoes, V_hat, symmetry)

Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, a unit vector representing the freestream.
"""
influence_matrix(colpoints, normals, horseshoes :: AbstractVector{Horseshoe}, V_hat :: SVector{3, <: Real}, symmetry :: Bool) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(r_i, horsie_j, n_i, V_hat, symmetry) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]


"""
    boundary_condition(velocities, normals)

Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)
