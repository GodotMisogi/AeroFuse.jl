## References
#==========================================================================================#

abstract type AbstractReferences end

struct References{T} <: AbstractReferences
    speed    :: T
    area     :: T
    span     :: T
    chord    :: T
    density  :: T
    location :: SVector{3,T}
end

References(V, S, b, c, ρ, ref :: AbstractVector{T}) where T <: Real = References{T}(V, S, b, c, ρ, ref)

References(; speed = 1., density = 1., span = 1., chord = 1., area = 1., location = zeros(3)) = References(speed, area, span, chord, density, location)

force_coefficient(force, refs :: References) = force_coefficient(force, dynamic_pressure(refs.density, refs.speed), refs.area)
moment_coefficient(moment, refs :: References) = moment_coefficient(moment, dynamic_pressure(refs.density, refs.speed), refs.area, refs.span, refs.chord)

rate_coefficient(fs :: Freestream, refs :: References) = rate_coefficient(fs.omega, refs.speed, refs.span, refs.chord)

## System
#==========================================================================================#

struct VLMSystem{M,N,R,S,P <: AbstractFreestream, Q <: AbstractReferences}
    vortices          :: M
    circulations      :: N 
    influence_matrix  :: R
    boundary_vector   :: S
    freestream        :: P
    reference         :: Q
end

# Miscellaneous
rate_coefficient(system :: VLMSystem) = rate_coefficient(system.freestream, system.reference)

# Velocities

surface_velocities(system :: VLMSystem) = surface_velocities(system.vortices, system.vortices, system.circulations, system.reference.speed * body_frame_velocity(system.freestream), system.freestream.omega)

## Forces

surface_forces(system :: VLMSystem, ::Geometry) = surface_forces(system.vortices, system.circulations, system.reference.speed * body_frame_velocity(system.freestream), system.freestream.omega, system.reference.density)

surface_forces(system :: VLMSystem, ::Stability) = geometry_to_stability_axes.(surface_forces(system, Geometry()), system.freestream.alpha)

surface_forces(system :: VLMSystem, ::Wind) = geometry_to_wind_axes.(surface_forces(system, Geometry()), system.freestream.alpha, system.freestream.beta)

surface_forces(system; axes :: AbstractAxisSystem = Geometry()) = surface_forces(system, axes)

## Moments

surface_moments(system :: VLMSystem, ::Geometry) = surface_moments(system.vortices, surface_forces(system, Geometry()), system.reference.location)

surface_moments(system :: VLMSystem, ::Stability) = surface_moments(system.vortices, stability_flip.(surface_forces(system, Stability())), system.reference.location)

surface_moments(system :: VLMSystem, ::Wind) = surface_moments(system.vortices, stability_flip.(surface_forces(system, Wind())), system.reference.location)

surface_moments(system; axes :: AbstractAxisSystem = Geometry()) = surface_moments(system, axes)

## Dynamics

function surface_dynamics(system :: VLMSystem)
    # Compute nearfield forces and moments
    surf_forces  = surface_forces(system)
    surf_moments = surface_moments(system.vortices, surf_forces, system.reference.location)

    surf_forces, surf_moments
end

surface_dynamics(system :: VLMSystem, ::Geometry) = surface_dynamics(system)

# TODO: Body axes
function surface_dynamics(system :: VLMSystem, ::Body)
end

function surface_dynamics(system :: VLMSystem, ::Stability)
    # Get angle of attack
    α = system.freestream.alpha

    # Compute surface forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to stability axes
    stability_forces  = geometry_to_stability_axes.(surface_forces, α)
    stability_moments = geometry_to_stability_axes.(stability_flip.(surface_moments), α)

    stability_forces, stability_moments
end

function surface_dynamics(system :: VLMSystem, ::Wind)
    # Get angles of attack and sideslip
    α = system.freestream.alpha
    β = system.freestream.beta

    # Compute nearfield forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to wind axes
    wind_forces  = @. geometry_to_wind_axes(surface_forces, α, β)
    wind_moments = @. geometry_to_wind_axes(stability_flip(surface_moments), α, β)

    wind_forces, wind_moments
end

surface_dynamics(system; axes :: AbstractAxisSystem = Wind()) = surface_dynamics(system, axes)

function surface_coefficients(system :: VLMSystem; axes :: AbstractAxisSystem = Wind()) 
    # Compute surface forces in whichever axes
    forces, moments = surface_dynamics(system, axes)

    # Compute coefficients
    CFs = force_coefficient.(forces, Ref(system.reference))
    CMs = moment_coefficient.(moments, Ref(system.reference))

    CFs, CMs
end

function nearfield_coefficients(system :: VLMSystem) 
    CFs, CMs = surface_coefficients(system; axes = Wind())
    ComponentVector(NamedTuple{keys(CFs)}([sum(CFs[key]); sum(CMs[key])] for key in valkeys(CFs)))
end

nearfield(system :: VLMSystem) = mapreduce(sum, vcat, surface_coefficients(system; axes = Wind()))

# Compute farfield forces in Trefftz plane
function farfield_forces(system :: VLMSystem)
    hs = system.vortices 
    Γs = system.circulations
    α  = system.freestream.alpha
    β  = system.freestream.beta
    V  = system.reference.speed
    ρ  = system.reference.density
    
    ComponentVector(NamedTuple{keys(hs)}(farfield_forces(Γs[comp], hs[comp], V, α, β, ρ) for comp in valkeys(hs)))
end

farfield_coefficients(system :: VLMSystem) = force_coefficient.(farfield_forces(system), dynamic_pressure(system.reference.density, system.reference.speed), system.reference.area)

farfield(system :: VLMSystem)  = let ff = farfield_coefficients(system); vec(sum(reshape(ff, 3, length(ff) ÷ 3), dims = 2)) end

# Made out of annoyance and boredom
function Base.show(io :: IO, sys :: VLMSystem)     
    println(io, "VLMSystem -")
    println(io, "Elements: ", length(sys.vortices), " ", eltype(sys.vortices))
    print(io, "Freestream: ")
    for fname in fieldnames(typeof(sys.freestream))
        print(io, fname, " = ", getfield(sys.freestream, fname), ", ")
    end
    print(io, "\nReferences: ")
    for fname in fieldnames(typeof(sys.reference))
        print(io, fname, " = ", getfield(sys.reference, fname), ", ")
    end
end