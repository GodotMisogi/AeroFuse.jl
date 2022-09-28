## Nearfield dynamics
#==========================================================================================#

kutta_joukowsky(ρ, V, l, Γ) = ρ * V × l * Γ

surface_velocity(h, horseshoes, Γs, U, Ω) = induced_trailing_velocity(bound_leg_center(h), horseshoes, Γs, U, Ω)

surface_velocities(hs_comp, horseshoes, Γs, U, Ω) = map(h -> surface_velocity(h, horseshoes, Γs, U, Ω), hs_comp)

"""
    surface_forces(hs_comp, Γ_comp, horseshoes, Γs, U, Ω, ρ)
    surface_forces(horseshoes, Γs, U, Ω, ρ)

Compute the surface forces via the local Kutta-Jowkowsky theorem.
 
For a given array of vortex type `hs_comp` and their associated vortex strengths ``Γ_c`` for which to compute the forces, the arrays of horseshoes and vortex strengths ``Γ``s  of the entire aircraft, the freestream flow vector ``U``, rotation rates ``Ω``, and a density ``ρ``, the velocities are evaluated at the midpoint of the bound leg of each horseshoe, excluding the contribution of the bound leg vortex.

The second variant simply sets `hs_comp = horseshoes` and `Γ_comp = Γs`.
"""
surface_forces(hs_comp, Γ_comp, horseshoes, Γs, U, Ω, ρ) = @timeit "KJForce" map((h, Γ) -> kutta_joukowsky(ρ, surface_velocity(h, horseshoes, Γs, U, Ω), bound_leg_vector(h), Γ), hs_comp, Γ_comp)

surface_forces(horseshoes, Γs, U, Ω, ρ) = surface_forces(horseshoes, Γs, horseshoes, Γs, U, Ω, ρ)

nearfield_drag(force, U) = -dot(force, normalize(U))

horseshoe_moment(horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force

surface_moments(horseshoes, forces, r_ref) = map((hs, f) -> horseshoe_moment(hs, f, r_ref), horseshoes, forces)

## Finite-core versions
#==========================================================================================#

# induced_trailing_velocity(r, horseshoes, Γs, U_hat, ε) = sum(x -> trailing_velocity(r, x[1], x[2], U_hat, ε), zip(horseshoes, Γs))
# induced_trailing_velocity(r, horseshoes, Γs, U, Ω, ε)          = induced_trailing_velocity(r, horseshoes, Γs, -normalize(U), ε) - (U + Ω × r)
# surface_velocity(hs, Γs, horseshoes, U, Ω, ε)          = @views induced_trailing_velocity(bound_leg_center(hs), horseshoes, Γs, U, Ω, ε)
# surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε)   = map(hs -> surface_velocity(hs, Γs, horseshoes, U, Ω, ε), hs_comp)

# surface_forces(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, ε) = map((Γ, hs, v) -> kutta_joukowsky(ρ, Γ, v, bound_leg_vector(hs)), Γ_comp, horseshoes, surface_velocities(hs_comp, Γs, horseshoes, U, Ω, ε))

# surface_forces(Γs, horseshoes, U, Ω, ρ, ε) = surface_forces(Γs, horseshoes, Γs, horseshoes, U, Ω, ρ, ε)

## Compressible variants using Prandtl-Glauert transformation
#==========================================================================================#

# surface_velocity(h, horseshoes, Γs, U, Ω, β) = prandtl_glauert_inverse_scale_velocity(surface_velocity(h, horseshoes, Γs, U, Ω), β)

# surface_velocities(hs_comp, horseshoes, Γs, U, Ω, β) = map(h -> surface_velocity(h, horseshoes, Γs, U, Ω, β), hs_comp)

# surface_forces(hs_comp, Γ_comp, horseshoes, Γs, U, Ω, ρ, β) = map((h, Γ) -> kutta_joukowsky(ρ, surface_velocity(h, horseshoes, Γs, U, Ω, β), bound_leg_vector(h), Γ), hs_comp, Γ_comp)
# surface_forces(horseshoes, Γs, U, Ω, ρ, β) = surface_forces(horseshoes, Γs, horseshoes, Γs, U, Ω, ρ, β)