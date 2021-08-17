function solve_aerodynamics!(Γ, system :: VLMSystem, state :: VLMState)
    # Update state velocity
    update_velocity!(state)

    # Assemble matrix system
    # @timeit "Influence Matrix" compute_influence_matrix!(system, state.velocity)
    # @timeit "Boundary Condition" compute_boundary_condition!(system, state.velocity, state.omega)
    @timeit "Generating AIC and RHS" generate_system!(system, state.velocity, state.omega) # Pre-allocated version for efficiency

    # Update circulations of system and surfaces
    @timeit "System Circulations"  system.circulations = Γ
    @timeit "Surface Circulations" update_circulations!(system)

    # Compute forces
    surfs = surfaces(system)
    @timeit "Surface Forces"  map(surf -> compute_surface_forces!(surf, system, state.velocity, state.omega, state.rho_ref), surfs)
    @timeit "Surface Moments" map(surf -> compute_surface_moments!(surf, state.r_ref), surfs)
    @timeit "Farfield Forces" map(surf -> compute_farfield_forces!(surf, state.speed, state.alpha, state.beta, state.rho_ref), surfs)

    nothing
end

# influence_matrix!(A, horseshoes, collocation_points, normals, V_hat, finite_core = false) = @einsum A[i,j] = influence_coefficient(horseshoes[j], collocation_points[i], normals[i], V_hat, finite_core)

# boundary_condition!(b, collocation_points, normals, V, Ω) = @einsum b[i] = dot(V + Ω × collocation_points[i], normals[i])

vlm_residual(hs, Γs, r, n, U_hat, Ω_hat) = dot(sum(x -> velocity(r, x[1], x[2], -U_hat), zip(hs, Γs)) - (U_hat + Ω_hat × r), n)
aerodynamic_residual!(R_A, hs, rs, nms, Γs, U_hat, Ω_hat) = @einsum R_A[i] = vlm_residual(hs, Γs, rs[i], nms[i], U_hat, Ω_hat)


## Nearfield improvement tests

function kutta_joukowsky(Γ_focus, hs_focus, Γs, horseshoes, U, Ω, ρ)
    U_hat = -normalize(U)
    funky(x) = trailing_velocity(x[1][1], x[1][2], x[2][1], x[2][2], U_hat, U, Ω, ρ)
    sum(funky, Base.Iterators.product(zip(hs_focus, Γ_focus), zip(horseshoes, Γs)))
end

trailing_velocity(hs, Γ, hs_all, Γ_all, U_hat, U, Ω, ρ) = ρ * (trailing_velocity(horseshoe_point(hs), hs_all, Γ_all, U_hat) - (U + Ω × horseshoe_point(hs))) × bound_leg_center(hs) * Γ