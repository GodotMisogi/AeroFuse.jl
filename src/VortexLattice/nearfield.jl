#-----------------------------Nearfield evaluations-----------------------#

"""
Placeholder.
"""
bound_leg_velocity(r, line :: Line, Γ, ε = 1e-8) = let a = r .- line.r1, b = r .- line.r2; on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ) end

"""
Placeholder.
"""
trailing_legs_velocities(r, line :: Line, Γ; direction = SVector(1., 0., 0.)) = @timeit "Trailing Leg" trailing_legs_velocities(r .- line.r1, r .- line.r2, Γ, direction)

"""
Placeholder.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r, bound_leg(horseshoe), Γ, direction = V)

"""
	midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluates the induced velocity by the trailing legs at the midpoint of a given Horseshoe ``r``, by summing over the velocities of Horseshoes with vortex strengths ``\\Gamma``s, rotation rates ``\\Omega``, and a freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r :: SVector{3, <: Real}, Ω :: SVector{3, <: Real}, horseshoes :: AbstractVector{Horseshoe}, Γs :: AbstractVector{<: Real}, U) = sum(trailing_velocity.(Ref(r), horseshoes, Γs, Ref(-normalize(U)))) .- U .- Ω × r

"""
	nearfield_forces(Γs, horseshoes, U, Ω, ρ)

Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths ``\\Gammas``s, a freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
nearfield_forces(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, U, Ω, ρ) = @timeit "Summing Forces" ρ * Γs .* midpoint_velocity.(bound_leg_center.(horseshoes), Ref(Ω), Ref(horseshoes), Ref(Γs), Ref(U)) .× bound_leg_vector.(horseshoes)

"""
	nearfield_drag(force, freestream)

Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Freestream.
"""
nearfield_drag(force, freestream :: Freestream) = dot(force, velocity(freestream) / freestream.mag)

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: AbstractVector{Horseshoe}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) ∈ zip(forces, horseshoes) ]

"""
Computes the nearfield forces and associated moments.
"""
function nearfield_dynamics(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225)
    @timeit "Forces" geom_forces = nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ)
    @timeit "Moments" geom_moments = moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end
