## Nearfield dynamics
#==========================================================================================#
kutta_joukowsky(ρ, Γ, V, l) = ρ * Γ * V × l

"""
Placeholder.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r1(r, horseshoe), r2(r, horseshoe), Γ, V)

"""
	midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluates the induced velocity by the trailing legs at a given location ``r``, by summing over the velocities of Horseshoes with vortex strengths ``\\Gamma``s, rotation rates ``\\Omega``, and a freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r, Ω, horseshoes :: AbstractVector{<: Horseshoe}, Γs, U) = sum(x -> trailing_velocity(r, x[1], x[2], -normalize(U)), zip(horseshoes, Γs)) - U - Ω × r

"""
	nearfield_forces(Γs, horseshoes, U, Ω, ρ)

Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths ``\\Gamma``s, a freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
function nearfield_forces(Γs, horseshoes :: AbstractVector{<: Horseshoe}, U, Ω, ρ)
    ρ_ref, Ω_ref, hs_ref, Γ_ref, U_ref = Ref(ρ), Ref(Ω), Ref(horseshoes), Ref(Γs), Ref(U)
    @. kutta_joukowsky(ρ_ref, Γs, midpoint_velocity(bound_leg_center(horseshoes), Ω_ref, hs_ref, Γ_ref, U_ref), bound_leg_vector(horseshoes))
end

sym_nearfield_forces(Γs, horseshoes :: AbstractVector{<: Horseshoe}, U, Ω, ρ) = [ nearfield_forces(Γs, reflect_xz.(horseshoes), U, Ω, ρ)[end:-1:1]; nearfield_forces(Γs, horseshoes, U, Ω, ρ) ]

"""
	nearfield_drag(force, freestream)

Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Freestream.
"""
nearfield_drag(force, freestream :: Freestream) = dot(force, velocity(freestream) / freestream.mag)

"""
Placeholder.
"""
horseshoe_moment(horseshoe :: Horseshoe, force, r_ref) = (bound_leg_center(horseshoe) - r_ref) × force

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: AbstractVector{<: Horseshoe}, forces, r_ref) = horseshoe_moment.(horseshoes, forces, Ref(r_ref))

"""
Computes the nearfield forces and associated moments.
"""
function nearfield_dynamics(Γs, horseshoes :: AbstractVector{<: Horseshoe}, freestream :: Freestream, r_ref, ρ, symmetry :: Bool)
    geom_forces = symmetry ? sym_nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ) : nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ)
    geom_moments = symmetry ? moments([ reflect_xz.(horseshoes[end:-1:1]); horseshoes ], geom_forces, r_ref) : moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end