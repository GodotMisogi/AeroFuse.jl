# Matrix setup
#==========================================================================================#

"""
	influence_coefficient(r, horseshoe, normal, V_hat, symmetry)

Compute the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
influence_coefficient(r, horseshoe :: Horseshoe, normal, V_hat) = dot(velocity(r, horseshoe, 1., V_hat), normal)

"""
	influence_matrix(colpoints, normals, horseshoes, V_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, and a unit vector representing the freestream.
"""
influence_matrix(colpoints, normals, horseshoes :: Vector{<: Horseshoe}, V_hat) = [ influence_coefficient(r_i, horsie_j, n_i, V_hat) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]

"""
	boundary_condition(velocities, normals)

Compute the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)