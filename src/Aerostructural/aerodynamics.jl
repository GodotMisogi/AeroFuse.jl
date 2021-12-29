# Induced velocities
induced_velocity(r, horseshoes, Γs, U_hat) = @timeit "Induced Velocity" @views sum(x -> velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

velocity(hs, Γs, r, n, U_hat, Ω_hat) = let vel = zeros(eltype(r), 3); @timeit "Velocity" dot(induced_velocity!(vel, r, hs, Γs, -U_hat) - (U_hat + Ω_hat × r), n) end

# Residual computation
aerodynamic_residual!(R_A, horseshoes, Γs, U_hat, Ω_hat) = @timeit "Aerodynamic Residual" R_A .= map(hs -> velocity(horseshoes, Γs, horseshoe_point(hs), horseshoe_normal(hs), U_hat, Ω_hat), horseshoes)

## USELESS? TESTING
#======================================#

function induced_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @views velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

# function surface_velocity!(vel, r, horseshoes, Γs, U, Ω)
#     U_hat = -normalize(U)
#     for i in eachindex(horseshoes)
#         vel += @views induced_trailing_velocity(r, horseshoes[i], Γs[i], U_hat) - (U + Ω × r)
#     end
    
#     vel
# end

# function evaluate_aerodynamics!(R, horseshoes, Γs, U, Ω, ρ, speed)
#     forces = Vector{typeof(U)}(undef, length(horseshoes))

#     U_hat = U / speed
#     Ω_hat = Ω / speed
#     for i in eachindex(horseshoes)
#         # Residual
#         R[i] = @views velocity(horseshoes, Γs / speed, horseshoe_point(horseshoes[i]), horseshoe_normal(horseshoes[i]), U_hat, Ω_hat)
        
#         # Forces
#         V_trail   = @views surface_velocity(horseshoes[i], Γs, horseshoes, U, Ω)
#         forces[i] = @views kutta_joukowsky(ρ, Γs[i], V_trail, bound_leg_vector(horseshoes[i]))
#     end

#     forces
# end


# function new_surface_forces(Γs_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, speed) 
#     U_hat  = -U / speed
#     forces = similar(Γs_comp, typeof(U))
#     for i in eachindex(hs_comp)
#         r_i = @views bound_leg_center(hs_comp[i])
#         l_i = @views bound_leg_vector(hs_comp[i])
#         vel = zeros(eltype(r_i), 3)
#         vel_fs = -(U + Ω × r_i)
#         for j in eachindex(horseshoes)
#             vel += @views trailing_velocity(r_i, horseshoes[j], Γs[j], U_hat) + vel_fs
#         end
#         forces[i] = @views kutta_joukowsky(ρ, Γs_comp[i], vel, l_i)
#    end

#    forces
# end
