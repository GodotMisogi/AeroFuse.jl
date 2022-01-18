## References
#==========================================================================================#

abstract type AbstractReferences end

"""
    References(V, ρ, μ, S, b, c, r)
    References(; speed, density, viscosity,
                 sound_speed, area, span, 
                 chord, location)
Define reference values with speed ``V``, density ``ρ``, _dynamic_ viscosity ``μ``, area ``S``, span ``b``, chord ``c``, location ``r`` for a vortex lattice analysis.

# Arguments
- `speed       :: Real         = 1.`: Speed (m/s)
- `density     :: Real         = 1.225`: Density (m)
- `viscosity   :: Real         = 1.5e-5`: Dynamic viscosity (m)
- `sound_speed :: Real         = 330.`: Speed of sound (m/s)
- `span        :: Real         = 1.`: Span length (m)
- `area        :: Real         = 1.`: Area (m²)
- `chord       :: Real         = 1.`: Chord length (m)
- `location    :: Vector{Real} = [0,0,0]`: Position
"""
struct References{T} <: AbstractReferences
    speed       :: T
    density     :: T
    viscosity   :: T
    sound_speed :: T
    area        :: T
    span        :: T
    chord       :: T
    location    :: SVector{3,T}
end

function References(V, ρ, μ, a, S, b, c, r)
    T = promote_type(eltype(V), eltype(ρ), eltype(μ), eltype(a), eltype(S), eltype(b), eltype(c), eltype(r))
    References{T}(V, ρ, μ, a, S, b, c, r) 
end

# References(; V, rho, mu, a, S, b, c, r) = References(V, rho, mu, a, S, b, c, r)

References(; speed = 1., density = 1.225, viscosity = 1.5e-5, sound_speed = 330., area = 1., span = 1., chord = 1., location = zeros(3)) = References(speed, density, viscosity, sound_speed, area, span, chord, location)

Base.broadcastable(refs :: References) = Ref(refs)

function Base.show(io :: IO, refs :: References)
    println(io, "\nReferences: ")
    for fname in fieldnames(typeof(refs))
        println(io, "    ", fname, " = ", getfield(refs, fname))
    end
end

kinematic_viscosity(refs :: References) = refs.viscosity / refs.density
mach_number(refs :: References) = refs.speed / refs.sound_speed

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
"""
    surface_velocities(system :: VLMSystem; 
                       axes   :: AbstractAxisSystem = Geometry())

Compute the velocities on the surface given the `VLMSystem` after performing an analysis, and the reference axis system.
"""
surface_velocities(system :: VLMSystem; axes = Geometry()) = surface_velocities(system, axes)

surface_velocities(system :: VLMSystem, ::Geometry) = surface_velocities(system.vortices, system.vortices, system.circulations, system.reference.speed * body_frame_velocity(system.freestream), system.freestream.omega)

surface_velocities(system :: VLMSystem, ::Body) = geometry_to_body_axes.(surface_velocities(system, Geometry()), system.freestream.alpha, system.freestream.beta)

surface_velocities(system :: VLMSystem, ::Stability) = geometry_to_stability_axes.(surface_velocities(system, Geometry()), system.freestream.alpha)

surface_velocities(system :: VLMSystem, ::Wind) = geometry_to_wind_axes.(surface_velocities(system, Geometry()), system.freestream.alpha, system.freestream.beta)

## Forces
"""
    surface_forces(system :: VLMSystem; 
                   axes   :: AbstractAxisSystem = Geometry())

Compute the forces on the surface given the `VLMSystem` after performing an analysis, and the reference axis system.
"""
surface_forces(system; axes :: AbstractAxisSystem = Geometry()) = surface_forces(system, axes)

surface_forces(system :: VLMSystem, ::Geometry) = surface_forces(system.vortices, system.circulations, system.reference.speed * body_frame_velocity(system.freestream), system.freestream.omega, system.reference.density)

surface_forces(system :: VLMSystem, ::Body) = geometry_to_body_axes.(surface_forces(system, Geometry()), system.freestream.alpha, system.freestream.beta)

surface_forces(system :: VLMSystem, ::Stability) = geometry_to_stability_axes.(surface_forces(system, Geometry()), system.freestream.alpha)

surface_forces(system :: VLMSystem, ::Wind) = geometry_to_wind_axes.(surface_forces(system, Geometry()), system.freestream.alpha, system.freestream.beta)

## Moments
"""
    surface_moments(system :: VLMSystem; 
                    axes   :: AbstractAxisSystem = Geometry())

Compute the moments on the surface given the `VLMSystem` after performing an analysis, and the reference axis system.
"""
surface_moments(system; axes :: AbstractAxisSystem = Geometry()) = surface_moments(system, axes)

surface_moments(system :: VLMSystem, ::Geometry) = surface_moments(system.vortices, surface_forces(system, Geometry()), system.reference.location)

surface_moments(system :: VLMSystem, ::Body) = surface_moments(system.vortices, surface_forces(system, Body()), system.reference.location)

surface_moments(system :: VLMSystem, ::Stability) = surface_moments(system.vortices, flip_xz.(surface_forces(system, Stability())), system.reference.location)

surface_moments(system :: VLMSystem, ::Wind) =  surface_moments(system.vortices, flip_xz.(surface_forces(system, Wind())), system.reference.location)


## Dynamics
"""
    surface_dynamics(system :: VLMSystem; 
                     axes   :: AbstractAxisSystem = Geometry())

Compute the forces and moments on the surface given the `VLMSystem` after performing an analysis, and the reference axis system.
"""
function surface_dynamics(system :: VLMSystem)
    # Compute nearfield forces and moments
    surf_forces  = surface_forces(system)
    surf_moments = surface_moments(system.vortices, surf_forces, system.reference.location)

    surf_forces, surf_moments
end

surface_dynamics(system :: VLMSystem, ::Geometry) = surface_dynamics(system)

function surface_dynamics(system :: VLMSystem, ::Body)
    # Compute surface forces and moments in geometry axes
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to body axes
    stability_forces  = @. geometry_to_body_axes(surface_forces)
    stability_moments = @. geometry_to_body_axes(surface_moments)

    stability_forces, stability_moments
end

function surface_dynamics(system :: VLMSystem, ::Stability)
    # Get angle of attack
    α = system.freestream.alpha

    # Compute surface forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to stability axes
    stability_forces  = @. geometry_to_stability_axes(surface_forces, α)
    stability_moments = @. geometry_to_stability_axes(flip_xz(surface_moments), α)

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
    wind_moments = @. geometry_to_wind_axes(flip_xz(surface_moments), α, β)

    wind_forces, wind_moments
end

surface_dynamics(system; axes :: AbstractAxisSystem = Wind()) = surface_dynamics(system, axes)

function surface_coefficients(system :: VLMSystem; axes :: AbstractAxisSystem = Wind()) 
    # Compute surface forces in whichever axes
    forces, moments = surface_dynamics(system, axes)

    # Compute coefficients
    CFs = @. force_coefficient(forces, system.reference)
    CMs = @. moment_coefficient(moments, system.reference)

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
    println(io, "Elements: ", length(sys.vortices), " ", eltype(sys.vortices), "\n")
    show(io, sys.freestream)
    show(io, sys.reference)
end