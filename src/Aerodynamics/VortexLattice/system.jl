## Cases
#==========================================================================================#

"""
    evaluate_case(horseshoe_panels :: Matrix{Panel3D}, normals, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), freestream speed `U`, angles of attack `α` and sideslip `β`,  reference density ``\\rho``, reference point ``r_\\text{ref}`` for moments, and reference values for area, chord, and span lengths.
"""
function evaluate_case(horseshoes :: Matrix{Horseshoe{T}}, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref) where T <: Real
    V = norm(U)

    # Solve system
    Γs, AIC, boco = solve_system(horseshoes[:], U, Ω)
    Γs = reshape(Γs, size(horseshoes))

    # Compute forces and moments
    surf_forces  = surface_forces(Γs, horseshoes, U, Ω, rho_ref)
    surf_moments = surface_moments(horseshoes, surf_forces, r_ref)
    far_force    = trefftz_forces(Γs, horseshoes, V, α, β, rho_ref)

    # Compute dynamic pressure
    q = dynamic_pressure(rho_ref, V)

    # Non-dimensional panel coefficients
    CFs = force_coefficient.(surf_forces, q, area_ref)
    CMs = moment_coefficient.(surf_moments, q, area_ref, span_ref, chord_ref)

    # Compute summed coefficients
    force, moment = sum(surf_forces), sum(surf_moments)

    # Transform near-field dynamics to wind axes
    trans_force  = body_to_wind_axes(force, α, β)
    trans_force  = [ nearfield_drag(force, U); trans_force[2:end] ]
    trans_moment = body_to_wind_axes(stability_flip(moment), α, β)

    # Compute coefficients
    nearfield_coeffs = aerodynamic_coefficients(trans_force, trans_moment, V, area_ref, span_ref, chord_ref, rho_ref)
    farfield_coeffs  = force_coefficient(far_force, q, area_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

abstract type AbstractReferences end

struct References{T} <: AbstractReferences
    area     :: T
    span     :: T
    chord    :: T
    density  :: T
    location :: SVector{3,T}
end

References(; density = 1., span = 1., chord = 1., area = 1., location :: AbstractVector{T} = zeros(3)) where T <: Real = References{T}(area, span, chord, density, location)

struct VLMSystem{M,N,R,S,P <: AbstractFreestream, Q <: AbstractReferences}
    horseshoes          :: M
    circulations        :: N 
    influence_matrix    :: R
    boundary_vector     :: S
    freestream          :: P
    reference           :: Q
end

abstract type AircraftAxes end

struct Global    <: AircraftAxes end
struct Body      <: AircraftAxes end
struct Stability <: AircraftAxes end
struct Wind      <: AircraftAxes end

surface_velocities(system :: VLMSystem) = surface_velocities(system.horseshoes, system.circulations, system.horseshoes, aircraft_velocity(system.freestream), system.freestream.omega)

surface_forces(system :: VLMSystem, ::Body) = surface_forces(system.circulations, system.horseshoes, aircraft_velocity(system.freestream), system.freestream.omega, system.reference.density)

surface_forces(system :: VLMSystem, ::Stability) = body_to_stability_axes.(surface_forces(system, Body()), system.freestream.alpha)

surface_forces(system :: VLMSystem, ::Wind) = body_to_wind_axes.(surface_forces(system, Body()), system.freestream.alpha, system.freestream.beta)

surface_forces(system; axes :: AircraftAxes = Body()) = surface_forces(system, axes)

surface_moments(system :: VLMSystem, ::Body) = surface_moments(system.horseshoes, surface_forces(system, Body()), system.reference.location)

surface_moments(system :: VLMSystem, ::Stability) = surface_moments(system.horseshoes, stability_flip.(surface_forces(system, Stability())), system.reference.location)

surface_moments(system :: VLMSystem, ::Wind) = surface_moments(system.horseshoes, stability_flip.(surface_forces(system, Wind())), system.reference.location)

surface_moments(system; axes :: AircraftAxes = Body()) = surface_moments(system, axes)

function surface_dynamics(system :: VLMSystem)
    hs = system.horseshoes
    r  = system.reference.location

    # Compute nearfield forces and moments
    surf_forces  = surface_forces(system)
    surf_moments = surface_moments(hs, surf_forces, r)

    surf_forces, surf_moments
end

# Body axes
surface_dynamics(system :: VLMSystem, ::Body) = surface_dynamics(system)

# Stability axes
function surface_dynamics(system :: VLMSystem, ::Stability)
    # Get angle of attack
    α = system.freestream.alpha

    # Compute surface forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to stability axes
    stability_forces  = body_to_stability_axes.(surface_forces, α)
    stability_moments = body_to_stability_axes.(stability_flip.(surface_moments), α)

    stability_forces, stability_moments
end

# Wind axes
function surface_dynamics(system :: VLMSystem, ::Wind)
    # Get angles of attack and sideslip
    α = system.freestream.alpha
    β = system.freestream.beta

    # Compute nearfield forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to wind axes
    wind_forces  = body_to_wind_axes.(surface_forces, α, β)
    wind_moments = body_to_wind_axes.(stability_flip.(surface_moments), α, β)

    wind_forces, wind_moments
end

surface_dynamics(system; axes :: AircraftAxes = Wind()) = surface_dynamics(system, axes)

function surface_coefficients(system :: VLMSystem; axes :: AircraftAxes = Wind()) 
    V = system.freestream.speed
    ρ = system.reference.density
    S = system.reference.area
    b = system.reference.span
    c = system.reference.chord

    # Compute surface forces in whichever axes
    forces, moments = surface_dynamics(system, axes)

    # Compute coefficients
    q   = dynamic_pressure(ρ, V)
    CFs = force_coefficient.(forces, q, S)
    CMs = moment_coefficient.(moments, q, S, b, c)

    CFs, CMs
end

# Compute farfield forces in Trefftz plane
function farfield_forces(system :: VLMSystem)
    hs = system.horseshoes 
    Γs = system.circulations
    V  = system.freestream.speed
    α  = system.freestream.alpha
    β  = system.freestream.beta
    ρ  = system.reference.density
    
    ComponentVector(NamedTuple{keys(hs)}(trefftz_forces(Γs[comp], hs[comp], V, α, β, ρ) for comp in valkeys(hs)))
end

function nearfield_coefficients(system :: VLMSystem) 
    CFs, CMs = surface_coefficients(system; axes = Wind())
    ComponentVector(NamedTuple{keys(CFs)}([sum(CFs[key]); sum(CMs[key])] for key in valkeys(CFs)))
end

farfield_coefficients(system :: VLMSystem) = force_coefficient.(farfield_forces(system), dynamic_pressure(system.reference.density, system.freestream.speed), system.reference.area)

nearfield(system :: VLMSystem) = mapreduce(sum, vcat, surface_coefficients(system; axes = Wind()))
farfield(system :: VLMSystem)  = let ff = farfield_coefficients(system); vec(sum(reshape(ff, 3, length(ff) ÷ 3), dims = 2)) end