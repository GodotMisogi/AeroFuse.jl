# Induced velocities
induced_velocity(r, horseshoes, Γs, U_hat) = @timeit "Induced Velocity" sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

velocity(hs, Γs, r, n, U_hat, Ω_hat) = @timeit "Velocity" dot(induced_velocity(r, hs, Γs, -U_hat) - (U_hat + Ω_hat × r), n)

# Residual computation
aerodynamic_residual!(R_A, horseshoes, Γs, U_hat, Ω_hat) = @timeit "Aerodynamic Residual" R_A .= map(hs -> velocity(horseshoes, Γs, horseshoe_point(hs), horseshoe_normal(hs), U_hat, Ω_hat), horseshoes)

function evaluate_aerodynamics!(R_A, Γ, all_horsies, U, Ω, ρ, speed) 
    @timeit "Surface Forces" all_forces = surface_forces(Γ, all_horsies, U, Ω, ρ)

    # new_forces = @views reshape(all_forces[1:length(new_horsies)], size(new_horsies))
    
    # Aerodynamic residuals
    @timeit "Aerodynamic Residual" aerodynamic_residual!(R_A, all_horsies, Γ / speed, U / speed, Ω / speed)

    all_forces
end