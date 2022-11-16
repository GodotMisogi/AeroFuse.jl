# Matrix setup
#==========================================================================================#

influence_coefficient(vor_i :: AbstractVortex, vor_j :: AbstractVortex) = dot(velocity(control_point(vor_j), vor_i, 1.), normal_vector(vor_j))

influence_coefficient(vor_i :: AbstractVortex, vor_j :: AbstractVortex, V_hat) = dot(velocity(control_point(vor_j), vor_i, 1., V_hat), normal_vector(vor_j))


"""
    influence_matrix(horseshoes)
    influence_matrix(horseshoes, u_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given an array of `Horseshoes` and the freestream direction ``̂u``.
"""
influence_matrix(horseshoes) = [ @timeit "Influence Coefficient" influence_coefficient(horsie_j, horsie_i) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]
influence_matrix(horseshoes, V_hat) = [ @timeit "Influence Coefficient" influence_coefficient(horsie_j, horsie_i, V_hat) for horsie_i ∈ horseshoes, horsie_j ∈ horseshoes ]

function influence_matrix!(A, horseshoes)
    for (i,j) in axes(A)
        A[i,j] = influence_coefficient(horseshoes[j], horseshoes[i])
    end
    return A
end

"""
    boundary_condition(horseshoes, U, Ω)

Assemble the boundary condition vector given an array of `Horseshoes`, the freestream velocity ``U``, and a quasi-steady rotation vector  ``Ω``.
"""
boundary_condition(horseshoes, U, Ω) = map(hs -> dot(U + Ω × control_point(hs), normal_vector(hs)), horseshoes)

# Added velocity for slipstream model
boundary_condition(horseshoes, U, Ups, Ω) = map((hs, Up) -> dot(U + Ω × control_point(hs) + Up, normal_vector(hs)), horseshoes, Ups)

# Matrix-free setup for nonlinear analyses
#==========================================================================================#

@views function induced_velocity(r, horseshoes, Γs, U_hat) 
    @timeit "Sum AD" vel = sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))
    
    # vel = zeros(eltype(Γs), 3)
    # vel1 = zeros(eltype(vel), 3)
    # for i in eachindex(horseshoes)
    #     velocity!(vel1, r, horseshoes[i], Γs[i], U_hat)
    #     vel += vel1
    # end

    return vel
end

@views function induced_trailing_velocity(r, horseshoes, Γs, U_hat) 
    @timeit "Sum AD" vel = sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

    # vel = zero(r)
    # @timeit "Sum AD" for i in eachindex(horseshoes)
    #     @timeit "Trailing Velocity" vel += trailing_velocity(r, horseshoes[i], Γs[i], U_hat)
    # end

    return vel
end

# In-place versions
@views function induced_velocity!(vel, r, horseshoes, Γs, U_hat)
    @timeit "Summing" for i in eachindex(horseshoes)
        vel += velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

# @views function induced_trailing_velocity!(vel, vel1, r, horseshoes, Γs, U_hat)
#     # for i in eachindex(horseshoes)
#         velocity!(vel1, r, horseshoes[i], Γs[i], U_hat)
#         vel += vel1
#     end

#     vel
# end

@views function induced_trailing_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += trailing_velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

function induced_velocity(r, hs, Γs, U, Ω) 
    # vel = zeros(eltype(Γs), 3)
    # vel1 = zeros(eltype(vel), 3)
    # induced_velocity!(vel, vel1, r, hs, Γs, -normalize(U)) - (U + Ω × r)
    vel = zero(r)
    @timeit "Induced Velocity" induced_velocity!(vel, r, hs, Γs, -normalize(U)) - (U + Ω × r)
    # @timeit "Induced Velocity" induced_velocity(r, hs, Γs, -normalize(U)) - (U + Ω × r)
end

function induced_trailing_velocity(r, horseshoes, Γs, U, Ω) 
    vel = zero(r)
    induced_trailing_velocity!(vel, r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)
    # induced_trailing_velocity(r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)
end

# Residual computations
residual(r, n, hs, Γs, U, Ω) = @timeit "Residual" dot(induced_velocity(r, hs, Γs, U, Ω), n)

solve_nonlinear(horseshoes, Γs, U_hat, Ω_hat) = map(hs -> residual(control_point(hs), normal_vector(hs), horseshoes, Γs, U_hat, Ω_hat), horseshoes)

solve_nonlinear!(R, horseshoes, Γs, U_hat, Ω_hat) = @timeit "Mapping" map!(hs -> residual(control_point(hs), normal_vector(hs), horseshoes, Γs, U_hat, Ω_hat), R, horseshoes)