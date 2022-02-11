# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(r, horseshoe, normal, u_hat, symmetry)

Compute the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``û`` at a point ``r`` projected to a normal vector.
"""
influence_coefficient(horseshoe :: Horseshoe, horseshoe_j, V_hat) = dot(velocity(horseshoe_point(horseshoe_j), horseshoe, 1., V_hat), horseshoe_normal(horseshoe_j))

"""
    influence_matrix(horseshoes, u_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given an array of `Horseshoes` and the freestream direction ``̂u``.
"""
influence_matrix(horseshoes, V_hat) = [ influence_coefficient(horsie_j, horsie_i, V_hat) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]

"""
    boundary_condition(horseshoes, U, Ω)

Assemble the boundary condition vector given an array of `Horseshoes`, the freestream velocity ``U``, and a quasi-steady rotation vector  ``Ω``.
"""
boundary_condition(horseshoes, U, Ω) = map(hs -> dot(U + Ω × horseshoe_point(hs), horseshoe_normal(hs)), horseshoes)

# Matrix-free setup for nonlinear analyses
#==========================================================================================#

induced_velocity(r, horseshoes, Γs, U_hat) = @views sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

# In-place version
function induced_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @views velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

induced_trailing_velocity(r, horseshoes, Γs, U_hat) = @views sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

# In-place version
function induced_trailing_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @views trailing_velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

function induced_velocity(r, hs, Γs, U, Ω) 
    # vel = zeros(eltype(r), 3)
    induced_velocity(r, hs, Γs, -normalize(U)) - (U + Ω × r)
end

function induced_trailing_velocity(r, horseshoes, Γs, U, Ω) 
    # vel = zeros(eltype(r), 3)
    induced_trailing_velocity(r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)
end

# Residual computations
function residual(r, n, hs, Γs, U, Ω) 
   dot(induced_velocity(r, hs, Γs, U, Ω), n)
end 

aerodynamic_residuals(horseshoes, Γs, U_hat, Ω_hat) = map(hs -> residual(horseshoe_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), horseshoes)

aerodynamic_residuals!(R, horseshoes, Γs, U_hat, Ω_hat) = map!(hs -> residual(horseshoe_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), R, horseshoes)