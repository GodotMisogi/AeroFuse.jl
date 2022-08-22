# Matrix setup
#==========================================================================================#

"""
    influence_coefficient(horseshoe_1 :: Horseshoe, horseshoe_2 :: Horseshoe)

Compute the influence coefficient of the first `Horseshoe` at the collocation point of the second `Horseshoe`.
"""
influence_coefficient(horseshoe :: Horseshoe, horseshoe_j) = dot(velocity(control_point(horseshoe_j), horseshoe, 1.), horseshoe_normal(horseshoe_j))
influence_coefficient(horseshoe :: Horseshoe, horseshoe_j, V_hat) = dot(velocity(control_point(horseshoe_j), horseshoe, 1., V_hat), horseshoe_normal(horseshoe_j))

"""
    influence_matrix(horseshoes)
    influence_matrix(horseshoes, u_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given an array of `Horseshoes` and the freestream direction ``̂u``.
"""
influence_matrix(horseshoes) = [ influence_coefficient(horsie_j, horsie_i) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]
influence_matrix(horseshoes, V_hat) = [ influence_coefficient(horsie_j, horsie_i, V_hat) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]

"""
    boundary_condition(horseshoes, U, Ω)

Assemble the boundary condition vector given an array of `Horseshoes`, the freestream velocity ``U``, and a quasi-steady rotation vector  ``Ω``.
"""
boundary_condition(horseshoes, U, Ω) = map(hs -> dot(U + Ω × control_point(hs), horseshoe_normal(hs)), horseshoes)

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

solve_nonlinear(horseshoes, Γs, U_hat, Ω_hat) = map(hs -> residual(control_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), horseshoes)

solve_nonlinear!(R, horseshoes, Γs, U_hat, Ω_hat) = map!(hs -> residual(control_point(hs), horseshoe_normal(hs), horseshoes, Γs, U_hat, Ω_hat), R, horseshoes)