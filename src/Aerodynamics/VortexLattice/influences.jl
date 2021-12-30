# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(r, horseshoe, normal, V_hat, symmetry)

Compute the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector.
"""
influence_coefficient(horseshoe :: Horseshoe, horseshoe_j, V_hat, finite_core = false) = dot(velocity(horseshoe_point(horseshoe_j), horseshoe, 1., V_hat, finite_core), horseshoe_normal(horseshoe_j))

"""
    influence_matrix(horseshoes, collocation_points, normals, V_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, and a unit vector representing the freestream.
"""
influence_matrix(horseshoes, V_hat, finite_core = false) = [ influence_coefficient(horsie_j, horsie_i, V_hat, finite_core) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]

"""
    boundary_condition(velocities, normals)

Compute the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = dot.(velocities, normals)

# Matrix-free setup for nonlinear analyses
#==========================================================================================#

induced_velocity(r, horseshoes, Γs, U_hat) = @timeit "Induced Velocity" @views sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

# In-place version
function induced_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @timeit "Induced Velocity" @views velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

# In-place version
function induced_trailing_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @timeit "Induced Trailing Velocity" @views trailing_velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

induced_trailing_velocity(r, horseshoes, Γs, U_hat) = @views sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

function induced_velocity(r, n, hs, Γs, U_hat, Ω_hat) 
    # vel = zeros(eltype(r), 3)
    @timeit "Velocity" dot(induced_velocity(r, hs, Γs, -U_hat) - (U_hat + Ω_hat × r), n) 
end

function midpoint_velocity(r, horseshoes, Γs, U, Ω) 
    # vel = zeros(eltype(r), 3)
    @timeit "Trailing Velocity" induced_trailing_velocity(r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)
end

# Residual computations
aerodynamic_residual(horseshoes, Γs, U_hat, Ω_hat) = map(hs -> induced_velocity(horseshoe_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), horseshoes)

aerodynamic_residual!(R, horseshoes, Γs, U_hat, Ω_hat) = @timeit "Aerodynamic Residual" map!(hs -> induced_velocity(horseshoe_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), R, horseshoes)

## Pre-allocated versions
#====================================================#

function matrix_assembly!(AIC, RHS, horseshoes, collocation_points, normals, V, Ω, finite_core = false)
    for inds ∈ CartesianIndices(AIC)
        i, j = inds.I
        AIC[inds] = dot(velocity(collocation_points[i], horseshoes[j], 1., -normalize(V), finite_core), normals[i])
        RHS[i]    = dot(V + Ω × collocation_points[i], normals[i])
    end
    nothing
end
    