## Nearfield dynamics
#==========================================================================================#

kutta_joukowsky(ρ, Γ, V, l) = ρ * Γ * V × l

"""

Evaluate the induced velocity at a given location ``r``, by summing over the trailing legs' velocities of Horseshoes with vortex strengths ``\\Gamma``s pointing in the direction ``\\hat U``.
"""
induced_trailing_velocity(r, horseshoes, Γs, U_hat) = @views sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

function induced_trailing_velocity!(vel, r, horseshoes, Γs, U_hat)
    for i in eachindex(horseshoes)
        vel += @views trailing_velocity(r, horseshoes[i], Γs[i], U_hat)
    end

    vel
end

"""
    midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluate the total velocity at a given location ``r`` by summing over the velocities induced by the trailing legs of Horseshoes with vortex strengths ``\\Gamma``s, of the rotation rates ``\\Omega``, and of the freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r, horseshoes, Γs, U, Ω) = induced_trailing_velocity(r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)

surface_velocity(hs, Γs, horseshoes, U, Ω) = midpoint_velocity(bound_leg_center(hs), horseshoes, Γs, U, Ω)

surface_velocities(hs_comp, Γs, horseshoes, U, Ω) = map(hs -> surface_velocity(hs, Γs, horseshoes, U, Ω), hs_comp)

"""
    surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ)

Compute the surface forces via the local Kutta-Jowkowski theorem given an array of horseshoes `hs_comp` to compute the forces on a component, their associated vortex strengths `Γ_comp`, the arrays of horseshoes and vortex strengths `Γs`  of the entire aircraft, the freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe, excluding the contribution of the bound leg.
"""

surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ) = map((Γ, hs) -> kutta_joukowsky(ρ, Γ, surface_velocity(hs, Γs, horseshoes, U, Ω), bound_leg_vector(hs)), Γ_comp, hs_comp)

surface_forces(Γs, horseshoes, U, Ω, ρ) = surface_forces(Γs, horseshoes, Γs, horseshoes, U, Ω, ρ)

nearfield_drag(force, U) = -dot(force, normalize(U))

horseshoe_moment(horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force

surface_moments(horseshoes, forces, r_ref) = map((hs, f) -> horseshoe_moment(hs, f, r_ref), horseshoes, forces)


## Finite-core versions
#==========================================================================================#

induced_trailing_velocity(r, horseshoes, Γs, U_hat, ε) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat, ε), zip(horseshoes, Γs))
midpoint_velocity(r, horseshoes, Γs, U, Ω, ε) = induced_trailing_velocity(r, horseshoes, Γs, -normalize(U), ε) - (U + Ω × r)
surface_velocity(hs, Γs, horseshoes, U, Ω, ε) = @views midpoint_velocity(bound_leg_center(hs), horseshoes, Γs, U, Ω, ε)
surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε) = map(hs -> surface_velocity(hs, Γs, horseshoes, U, Ω, ε), hs_comp)

surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, ε) = map((Γ, hs, v) -> kutta_joukowsky(ρ, Γ, v, bound_leg_vector(hs)), Γ_comp, horseshoes, surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε))

surface_forces(Γs, horseshoes, U, Ω, ρ, ε) = surface_forces(Γs, horseshoes, Γs, horseshoes, U, Ω, ρ, ε)