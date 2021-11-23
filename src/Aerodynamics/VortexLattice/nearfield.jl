## Nearfield dynamics
#==========================================================================================#

kutta_joukowsky(ρ, Γ, V, l) = ρ * Γ * V × l

"""

Evaluate the induced velocity at a given location ``r``, by summing over the trailing legs' velocities of Horseshoes with vortex strengths ``\\Gamma``s pointing in the direction ``\\hat U``.
"""
induced_trailing_velocity(r, horseshoes, Γs, U_hat) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat), zip(horseshoes, Γs))

"""
    midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluate the total velocity at a given location ``r`` by summing over the velocities induced by the trailing legs of Horseshoes with vortex strengths ``\\Gamma``s, of the rotation rates ``\\Omega``, and of the freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r, horseshoes, Γs, U, Ω) = induced_trailing_velocity(r, horseshoes, Γs, -normalize(U)) - (U + Ω × r)

panel_velocities(hs_comp, Γs, horseshoes, U, Ω) = map(x -> midpoint_velocity(bound_leg_center(x), horseshoes, Γs, U, Ω), hs_comp)

"""
    nearfield_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ)

Compute the nearfield forces via the local Kutta-Jowkowski theorem given an array of horseshoes `hs_comp` to compute the forces on a component, their associated vortex strengths `Γ_comp`, the arrays of horseshoes and vortex strengths `Γs`  of the entire aircraft, the freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe, excluding the contribution of the bound leg.
"""
nearfield_forces(Γs, horseshoes, U, Ω, ρ) = map((Γ, hs, v) -> kutta_joukowsky(ρ, Γ, v, bound_leg_vector(hs)), Γs, horseshoes, panel_velocities(horseshoes, Γs, horseshoes, U, Ω))

nearfield_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ) = map((Γ, hs, v) -> kutta_joukowsky(ρ, Γ, v, bound_leg_vector(hs)), Γ_comp, horseshoes, panel_velocities(hs_comp, Γs, horseshoes, U, Ω))

nearfield_drag(force, U) = -dot(force, normalize(U))

horseshoe_moment(horseshoe :: Horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force

nearfield_moments(horseshoes, forces, r_ref) = map((hs, f) -> horseshoe_moment(hs, f, r_ref), horseshoes, forces)