# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(r, horseshoe, panel_normal, V_hat, symmetry)

Computes the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
influence_coefficient(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, panel_normal :: SVector{3, <: Real}, V_hat :: SVector{3, <: Real}) = dot(velocity(r, horseshoe, 1., V_hat), panel_normal)

"""
    influence_matrix(colpoints, normals, horseshoes, V_hat, symmetry)

Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, a unit vector representing the freestream.
"""
influence_matrix(colpoints, normals, horseshoes :: AbstractVector{Horseshoe}, V_hat :: SVector{3, <: Real}) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(r_i, horsie_j, n_i, V_hat) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]


"""
    boundary_condition(velocities, normals)

Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)
