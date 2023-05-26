# Matrix setup
#==========================================================================================#

influence_coefficient(vor_i :: AbstractVortex, vor_j :: AbstractVortex) = dot(velocity(control_point(vor_j), vor_i, 1.), normal_vector(vor_j))

influence_coefficient(vor_i :: AbstractVortex, vor_j :: AbstractVortex, V_hat) = dot(velocity(control_point(vor_j), vor_i, 1., V_hat), normal_vector(vor_j))


"""
    influence_matrix(vortices)
    influence_matrix(vortices, u_hat)

Assemble the Aerodynamic Influence Coefficient (AIC) matrix given an array of `AbstractVortex` and the freestream direction ``û``.
"""
influence_matrix(vortices) = [ influence_coefficient(horsie_j, horsie_i) for horsie_i ∈ vortices, horsie_j ∈ vortices ]
influence_matrix(vortices, V_hat) = [ influence_coefficient(horsie_j, horsie_i, V_hat) for horsie_i ∈ vortices, horsie_j ∈ vortices ]

function influence_matrix!(A, vortices)
    for (i,j) in axes(A)
        A[i,j] = influence_coefficient(vortices[j], vortices[i])
    end
    return A
end

"""
    boundary_condition(vortices, U, Ω)

Assemble the boundary condition vector given an array of `AbstractVortex`, the freestream velocity ``U``, and a quasi-steady rotation vector  ``Ω``.
"""
boundary_condition(vortices, U, Ω) = map(hs -> dot(U + Ω × control_point(hs), normal_vector(hs)), vortices)

# Added velocity for slipstream model
boundary_condition(vortices, U, Ups, Ω) = map((hs, Up) -> dot(U + Ω × control_point(hs) + Up, normal_vector(hs)), vortices, Ups)

# Matrix-free setup for nonlinear analyses
#==========================================================================================#

induced_velocity(r, vortices, Γs, U_hat) = sum(x -> velocity(r, x[1], x[2], U_hat), zip(vortices, Γs))

induced_trailing_velocity(r, vortices, Γs, U_hat) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(vortices, Γs))

# In-place versions
@views function induced_velocity!(vel, r, vortices, Γs, U_hat)
    for i in eachindex(vortices)
        vel += velocity(r, vortices[i], Γs[i], U_hat)
    end

    return vel
end

@views function induced_trailing_velocity!(vel, r, vortices, Γs, U_hat)
    for i in eachindex(vortices)
        vel += trailing_velocity(r, vortices[i], Γs[i], U_hat)
    end

    return vel
end

function induced_velocity(r, hs, Γs, U, Ω)
    vel = zero(r)
    induced_velocity!(vel, r, hs, Γs, -normalize(U)) - (U + Ω × r)
    # @timeit "Induced Velocity" induced_velocity(r, hs, Γs, -normalize(U)) - (U + Ω × r)
end

function induced_trailing_velocity(r, vortices, Γs, U, Ω) 
    vel = zero(r)
    induced_trailing_velocity!(vel, r, vortices, Γs, -normalize(U)) - (U + Ω × r)
    # induced_trailing_velocity(r, vortices, Γs, -normalize(U)) - (U + Ω × r)
end

# Residual computations
residual(r, n, hs, Γs, U, Ω) = dot(induced_velocity(r, hs, Γs, U, Ω), n)

solve_nonlinear(vortices, Γs, U_hat, Ω_hat) = map(hs -> residual(control_point(hs), normal_vector(hs), vortices, Γs, U_hat, Ω_hat), vortices)

solve_nonlinear!(R, vortices, Γs, U_hat, Ω_hat) = map!(hs -> residual(control_point(hs), normal_vector(hs), vortices, Γs, U_hat, Ω_hat), R, vortices)