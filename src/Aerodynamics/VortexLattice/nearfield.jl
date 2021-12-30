## Nearfield dynamics
#==========================================================================================#

kutta_joukowsky(ρ, V, l, Γ) = ρ * V × l * Γ

surface_velocity(h, Γs, horseshoes, U_hat, Ω_hat) = midpoint_velocity(bound_leg_center(h), horseshoes, Γs, U_hat, Ω_hat)

surface_velocities(hs_comp, Γs, horseshoes, U_hat, Ω_hat) = map(h -> surface_velocity(h, Γs, horseshoes, U_hat, Ω_hat), hs_comp)

"""
    surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ)

Compute the surface forces via the local Kutta-Jowkowski theorem given an array of horseshoes `hs_comp` to compute the forces on a component, their associated vortex strengths `Γ_comp`, the arrays of horseshoes and vortex strengths `Γs`  of the entire aircraft, the freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe, excluding the contribution of the bound leg.
"""

surface_forces(hs_comp, Γ_comp, horseshoes, Γs, speed, U_hat, Ω_hat, ρ) = map((Γ, h) -> kutta_joukowsky(ρ, speed * surface_velocity(h, Γs, horseshoes, U_hat, Ω_hat), bound_leg_vector(h), speed * Γ), Γ_comp, hs_comp)

surface_forces(horseshoes, Γs, speed, U_hat, Ω_hat, ρ) = surface_forces(horseshoes, Γs, horseshoes, Γs, speed, U_hat, Ω_hat, ρ)

nearfield_drag(force, U) = -dot(force, normalize(U))

horseshoe_moment(horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force

surface_moments(horseshoes, forces, r_ref) = map((hs, f) -> horseshoe_moment(hs, f, r_ref), horseshoes, forces)


## Finite-core versions
#==========================================================================================#

induced_trailing_velocity(r, horseshoes, Γs, U_hat, ε) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat, ε), zip(horseshoes, Γs))
midpoint_velocity(r, horseshoes, Γs, U, Ω, ε)          = induced_trailing_velocity(r, horseshoes, Γs, -U, ε) - (U + Ω × r)
surface_velocity(hs, Γs, horseshoes, U, Ω, ε)          = @views midpoint_velocity(bound_leg_center(hs), horseshoes, Γs, U, Ω, ε)
surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε)   = map(hs -> surface_velocity(hs, Γs, horseshoes, U, Ω, ε), hs_comp)

surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, ε) = map((Γ, hs, v) -> kutta_joukowsky(ρ, Γ, v, bound_leg_vector(hs)), Γ_comp, horseshoes, surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε))

surface_forces(Γs, horseshoes, U, Ω, ρ, ε) = surface_forces(Γs, horseshoes, Γs, horseshoes, U, Ω, ρ, ε)