function solve_aerodynamics!(Γ, system :: VLMSystem, state :: VLMState)
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

# Induced velocities
induced_velocity(r, horseshoes, Γs, U_hat) = sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

velocity(hs, Γs, r, n, U_hat, Ω_hat) = dot(induced_velocity(r, hs, Γs, -U_hat) - (U_hat + Ω_hat × r), n)

# Residual computation
aerodynamic_residual!(R_A, hs, rs, nms, Γs, U_hat, Ω_hat) = R_A .= velocity.(Ref(hs), Ref(Γs), rs, nms, Ref(U_hat), Ref(Ω_hat))

# HEAVILY ALLOCATING AERODYNAMIC RESIDUAL COMPUTATIONS
# 
# Single surface:
# @timeit "Influence Matrix" inf_mat = influence_matrix(new_horsies[:], horseshoe_point.(new_horsies[:]), new_normies[:], -normalize(U), false)
# @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(new_horsies[:], U, Ω), new_normies[:])
# @timeit "Aerodynamic Residual" linear_residual!(R_A, inf_mat, Γ, boco)
#
# Multiple surfaces:
# @timeit "Normalizing Velocity" U_hat = normalize(U)
# @timeit "Influence Matrix" inf_mat = influence_matrix(all_horsies, horseshoe_point.(all_horsies), normies[:], -U_hat, false)
# @timeit "Boundary Condition" boco    = boundary_condition(quasi_steady_freestream(all_horsies, U_hat, Ω / speed), normies[:])
# @timeit "Aerodynamic Residual" linear_residual!(R_A, inf_mat, Γ / speed, boco)