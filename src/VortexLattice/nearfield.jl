## Nearfield dynamics
#==========================================================================================#
kutta_joukowsky(ρ, Γ, V, l) = ρ * Γ * V × l

"""
Placeholder.
"""
trailing_velocity(r :: SVector{3,<: Real}, horseshoe :: Horseshoe, Γ :: Real, V :: SVector{3,<: Real}) = 
    let line = bound_leg(horseshoe)
        @timeit "a" a = r - point1(line)
        @timeit "b" b = r - point2(line); 
        @timeit "Trailing Velocity" trailing_legs_velocities(a, b, Γ, V) end

"""
	midpoint_velocity(r, Ω, horseshoes, Γs, U)

Evaluates the induced velocity by the trailing legs at a given location ``r``, by summing over the velocities of Horseshoes with vortex strengths ``\\Gamma``s, rotation rates ``\\Omega``, and a freestream flow vector ``U`` in the aircraft reference frame.
"""
midpoint_velocity(r :: SVector{3, <: Real}, Ω :: SVector{3,<: Real}, horseshoes :: AbstractVector{Horseshoe}, Γs :: AbstractVector{<: Real}, U :: SVector{3,<: Real}) = @timeit "Midpoint Velocity" sum(trailing_velocity.(Ref(r), horseshoes, Γs, Ref(-normalize(U)))) - U - Ω × r

# kutta_joukowsky(r :: SVector{3, <: Real}, Γ :: Real, l :: SVector{3, <: Real}, Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, U :: SVector{3,<: Real}, Ω :: SVector{3,<: Real}, ρ :: Real) = kutta_joukowsky(ρ, Γ, midpoint_velocity(r, Ω, horseshoes, Γs, U), l)

"""
	nearfield_forces(Γs, horseshoes, U, Ω, ρ)

Computes the nearfield forces via the local Kutta-Jowkowski given an array of horseshoes, their associated vortex strengths ``\\Gamma``s, a freestream flow vector ``U``, rotation rates ``\\Omega``, and a density ``\\rho``. The velocities are evaluated at the midpoint of the bound leg of each horseshoe.
"""
nearfield_forces(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, U :: SVector{3,<: Real}, Ω :: SVector{3,<: Real}, ρ :: Real) = @timeit "Summing Forces" forces = kutta_joukowsky.(Ref(ρ), Γs, midpoint_velocity.(bound_leg_center.(horseshoes), Ref(Ω), Ref(horseshoes), Ref(Γs), Ref(U)), bound_leg_vector.(horseshoes))

sym_nearfield_forces(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, U :: SVector{3,<: Real}, Ω :: SVector{3,<: Real}, ρ :: Real) = [ nearfield_forces(Γs, reflect_xz.(horseshoes), U, Ω, ρ)[end:-1:1]; nearfield_forces(Γs, horseshoes, U, Ω, ρ) ]

"""
	nearfield_drag(force, freestream)

Computes the near-field drag given the sum of the local Kutta-Jowkowski forces and a Freestream.
"""
nearfield_drag(force :: SVector{3,<: Real}, freestream :: Freestream) = dot(force, velocity(freestream) / freestream.mag)

"""
Placeholder.
"""
horseshoe_moment(horseshoe :: Horseshoe, force :: SVector{3,<: Real}, r_ref :: SVector{3,<: Real}) = (bound_leg_center(horseshoe) - r_ref) × force

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: AbstractVector{Horseshoe}, forces, r_ref :: SVector{3,<: Real}) = horseshoe_moment.(horseshoes, forces, Ref(r_ref))

"""
Computes the nearfield forces and associated moments.
"""
function nearfield_dynamics(Γs :: AbstractVector{<: Real}, horseshoes :: AbstractVector{Horseshoe}, freestream :: Freestream, r_ref :: SVector{3,<: Real}, ρ :: Real, symmetry :: Bool)
    @timeit "Forces" geom_forces = symmetry ? sym_nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ) : nearfield_forces(Γs, horseshoes, aircraft_velocity(freestream), freestream.Ω, ρ)
    @timeit "Moments" geom_moments = symmetry ? moments([ reflect_xz.(horseshoes[end:-1:1]); horseshoes ], geom_forces, r_ref) : moments(horseshoes, geom_forces, r_ref)

    geom_forces, geom_moments
end