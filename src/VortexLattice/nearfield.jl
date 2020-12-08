#-----------------------------Nearfield evaluations-----------------------#

"""
Placeholder.
"""
bound_leg_velocity(r, line :: Line, Γ, ε = 1e-8) = let a = r .- line.r1, b = r .- line.r2; on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ) end

"""
Placeholder.
"""
trailing_legs_velocities(r, line :: Line, Γ; direction = SVector(1., 0., 0.)) = @timeit "Trailing Leg" let a = r .- line.r1, b = r .- line.r2; trailing_legs_velocities(a, b, Γ, direction) end

"""
Placeholder.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r, horseshoe.bound_leg, Γ, direction = V)

"""
Evaluates the induced velocity by the trailing legs at the midpoint of a given Horseshoe `r`, by summing over the velocities of Horseshoes with vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the aircraft reference frame.
"""
midpoint_velocity(r :: SVector{3, <: Real}, Ω :: SVector{3, <: Real}, horseshoes :: AbstractVector{Horseshoe}, Γs :: AbstractVector{<: Real}, U) = sum(trailing_velocity(r, horseshoe, Γ, -normalize(U)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .- U .- Ω × r

"""
Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths Γs, a Freestream, and a density ρ. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
function nearfield_forces(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, U, Ω, ρ)
    Γ_shoes = zip(Γs, horseshoes)
    @timeit "Summing Forces" [ 
      let r_i = bound_leg_center(horseshoe_i), l_i = bound_leg_vector(horseshoe_i); 
      ρ * Γ_i * midpoint_velocity(r_i, Ω, horseshoes, Γs, U) × l_i end 
      for (Γ_i, horseshoe_i) ∈ Γ_shoes 
    ]
end

"""
Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Freestream.
"""
nearfield_drag(force, V :: Freestream) = dot(force, velocity(V) / V.mag)

"""
Computes the nearfield forces and associated moments.
"""
function nearfield_dynamics(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225)
    @timeit "Forces" geom_forces = nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ)
    @timeit "Moments" geom_moments = moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end
