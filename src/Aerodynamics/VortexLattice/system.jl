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
- `area        :: Real         = 1.`: Area (m²)
- `span        :: Real         = 1.`: Span length (m)
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

dynamic_pressure(refs :: References) = 1/2 * refs.density * refs.speed^2

kinematic_viscosity(refs :: References) = refs.viscosity / refs.density
mach_number(refs :: References) = refs.speed / refs.sound_speed

force_coefficient(force, refs :: References) = force_coefficient(force, dynamic_pressure(refs), refs.area)
moment_coefficient(moment, refs :: References) = moment_coefficient(moment, dynamic_pressure(refs), refs.area, refs.span, refs.chord)

rate_coefficient(fs :: Freestream, refs :: References) = rate_coefficient(fs.omega, refs.speed, refs.span, refs.chord)

## System
#==========================================================================================#

abstract type AbstractVortexLatticeSystem end

"""
    VortexLatticeSystem

A system consisting of the relevant variables for a vortex lattice analysis for post-processing.

# Arguments
The accessible fields are:
- `vortices`:  The array of vortices (currently `Horseshoe`).
- `circulations`: The circulation strengths of the vortices obtained by solving the linear system.
- `influence_matrix`: The influence matrix of the linear system.
- `boundary_vector`: The boundary condition corresponding to the right-hand-side of the linear system.
- `freestream :: Freestream`: The freestream conditions.
- `reference :: References`: The reference values.
"""
struct VortexLatticeSystem{
    M <: DenseArray{<: AbstractVortex},
    N <: DenseArray{<: Real},
    R,
    S,
    P <: AbstractFreestream,
    Q <: AbstractReferences,
    T <: AbstractAxisSystem} <: AbstractVortexLatticeSystem
    vortices          :: M
    circulations      :: N 
    influence_matrix  :: R
    boundary_vector   :: S
    freestream        :: P
    reference         :: Q
    compressible      :: Bool
    axes              :: T
end

"""
    VortexLatticeSystem(
        aircraft, 
        fs :: Freestream, 
        refs :: References, 
        compressible = false, 
        axes = Geometry()
    )

Construct a `VortexLatticeSystem` for analyzing inviscid aerodynamics of an aircraft (must be a `ComponentArray` of `Horseshoe`s or `VortexRing`s) with `Freestream` conditions and `References` for non-dimensionalization. Options are provided for compressibility corrections via the Prandtl-Glauert transformation (false by default) and axis system for computing velocities and forces (`Geometry` by default).
"""
function VortexLatticeSystem(aircraft, fs :: Freestream, refs :: References, compressible = false, axes = Geometry())

    # Mach number bound checks
    M = mach_number(refs)
    @assert M < 1.  "Only compressible subsonic flow conditions (M < 1) are valid!"
    if M > 0.7 @warn "Results in transonic to sonic flow conditions (0.7 < M < 1) are most likely incorrect!" end
    if M > 0.3 && !compressible @warn "Compressible regime (M > 0.3) but compressibility correction is off, be wary of the analysis!" end

    # (Prandtl-Glauert ∘ Wind axis) transformation
    if compressible
        β_pg = √(1 - M^2)
        ac = @. prandtl_glauert_scale_coordinates(geometry_to_wind_axes(aircraft, fs), β_pg)
    else
        β_pg = 1
        ac = @. geometry_to_wind_axes(aircraft, fs)
    end

    # Quasi-steady freestream velocity
    U = geometry_to_wind_axes(velocity(fs, Body()), fs)
    Ω = geometry_to_wind_axes(fs.omega, fs) / refs.speed

    # Solve system
    Γs, AIC, boco = solve_linear(ac, U, Ω)

    return VortexLatticeSystem(aircraft, refs.speed * Γs / β_pg^2, AIC, boco, fs, refs, compressible, axes)
end

# Miscellaneous
rate_coefficient(system :: VortexLatticeSystem) = rate_coefficient(system.freestream, system.reference)

## THINK ABOUT USING ONLY WIND AXES FOR PG-TRANSFORMATION AND MAPPING BACK

## Velocities
"""
    surface_velocities(system :: VortexLatticeSystem; 
                       axes   :: AbstractAxisSystem = Geometry())

Compute the induced velocities for all components of the `VortexLatticeSystem` in the reference axis system.
"""
surface_velocities(system :: VortexLatticeSystem; axes = Geometry()) = surface_velocities(system, axes)

surface_velocities(system :: VortexLatticeSystem, ::Geometry) = surface_velocities(system.vortices, system.vortices, system.circulations, system.reference.speed * velocity(system.freestream, Body()), system.freestream.omega)

surface_velocities(system :: VortexLatticeSystem, ::Body) = geometry_to_body_axes.(surface_velocities(system, Geometry()), system.freestream.alpha, system.freestream.beta)

surface_velocities(system :: VortexLatticeSystem, ::Stability) = geometry_to_stability_axes.(surface_velocities(system, Geometry()), system.freestream.alpha)

surface_velocities(system :: VortexLatticeSystem, ::Wind) = geometry_to_wind_axes.(surface_velocities(system, Geometry()), system.freestream.alpha, system.freestream.beta)

## Forces
"""
    surface_forces(system :: VortexLatticeSystem; 
                   axes   :: AbstractAxisSystem = Geometry())

Compute the forces for all components of the `VortexLatticeSystem` in the reference axis system.
"""
surface_forces(system; axes :: AbstractAxisSystem = Geometry()) = surface_forces(system, axes)

surface_forces(system :: VortexLatticeSystem, ::Geometry) = surface_forces(system.vortices, system.circulations, system.reference.speed * velocity(system.freestream, Body()), system.freestream.omega, system.reference.density)

surface_forces(system :: VortexLatticeSystem, ::Body) = geometry_to_body_axes.(surface_forces(system, Geometry()), system.freestream.alpha, system.freestream.beta)

surface_forces(system :: VortexLatticeSystem, ::Stability) = geometry_to_stability_axes.(surface_forces(system, Geometry()), system.freestream.alpha)

surface_forces(system :: VortexLatticeSystem, ::Wind) = geometry_to_wind_axes.(surface_forces(system, Geometry()), system.freestream.alpha, system.freestream.beta)

## Moments
"""
    surface_moments(system :: VortexLatticeSystem; 
                    axes   :: AbstractAxisSystem = Geometry())

Compute the moments for all components of the `VortexLatticeSystem` in the reference axis system.
"""
surface_moments(system; axes :: AbstractAxisSystem = Geometry()) = surface_moments(system, axes)

surface_moments(system :: VortexLatticeSystem, ::Geometry) = surface_moments(system.vortices, surface_forces(system, Geometry()), system.reference.location)

surface_moments(system :: VortexLatticeSystem, ::Body) = surface_moments(system.vortices, surface_forces(system, Body()), geometry_to_body_axes(system.reference.location))

surface_moments(system :: VortexLatticeSystem, ::Stability) = surface_moments(system.vortices, flip_xz.(surface_forces(system, Stability())), geometry_to_stability_axes(system.reference.location, system.freestream.alpha))

surface_moments(system :: VortexLatticeSystem, ::Wind) =  surface_moments(system.vortices, flip_xz.(surface_forces(system, Wind())), geometry_to_wind_axes(system.reference.location, system.freestream.alpha, system.freestream.beta))


## Dynamics
function surface_dynamics(system :: VortexLatticeSystem)
    # Compute nearfield forces and moments
    surf_forces  = surface_forces(system)
    surf_moments = surface_moments(system.vortices, surf_forces, system.reference.location)

    surf_forces, surf_moments
end

"""
    surface_dynamics(system :: VortexLatticeSystem; 
                     axes   :: AbstractAxisSystem = Geometry())

Compute the forces and moments for all components of the `VortexLatticeSystem` in the reference axis system.
"""
surface_dynamics(system; axes :: AbstractAxisSystem = Wind()) = surface_dynamics(system, axes)

surface_dynamics(system :: VortexLatticeSystem, ::Geometry) = surface_dynamics(system)

function surface_dynamics(system :: VortexLatticeSystem, ::Body)
    # Compute surface forces and moments in geometry axes
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to body axes
    stability_forces  = @. geometry_to_body_axes(surface_forces)
    stability_moments = @. geometry_to_body_axes(surface_moments)

    stability_forces, stability_moments
end

function surface_dynamics(system :: VortexLatticeSystem, ::Stability)
    # Get angle of attack
    α = system.freestream.alpha

    # Compute surface forces and moments
    surface_forces, surface_moments = surface_dynamics(system)

    # Transform to stability axes
    stability_forces  = @. geometry_to_stability_axes(surface_forces, α)
    stability_moments = @. geometry_to_stability_axes(flip_xz(surface_moments), α)

    stability_forces, stability_moments
end

function surface_dynamics(system :: VortexLatticeSystem, ::Wind)
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

"""
    surface_coefficients(system :: VortexLatticeSystem; 
                         axes   :: AbstractAxisSystem = Wind())

Compute the force and moment coefficients on the surface given the `VortexLatticeSystem` in the reference axis system.
"""
function surface_coefficients(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Wind()) 
    # Compute surface forces in whichever axes
    forces, moments = surface_dynamics(system, axes)
    refs = system.reference

    # Compute coefficients
    CFs = @. force_coefficient(forces, refs)
    CMs = @. moment_coefficient(moments, dynamic_pressure(refs), refs.area, refs.span, refs.chord)

    CFs, CMs
end

const NF_COEFFS = @SLArray (6) (:CX,:CY,:CZ,:Cl,:Cm,:Cn)
const FF_COEFFS = @SLArray (3) (:CDi,:CY,:CL)

"""
    nearfield_coefficients(system :: VortexLatticeSystem)

Compute the force and moment coefficients in **wind axes** for all components of the `VortexLatticeSystem`.
"""
function nearfield_coefficients(system :: VortexLatticeSystem) 
    CFs, CMs = surface_coefficients(system; axes = Wind())
 
    @views NamedTuple(key => NF_COEFFS(sum(CFs[key])..., sum(CMs[key])...) for key in keys(CFs))
end

"""
    nearfield(system :: VortexLatticeSystem)

Compute the **total** force and moment coefficients in **wind axes** for all components of the `VortexLatticeSystem`.
"""
function nearfield(system :: VortexLatticeSystem)
    CX, CY, CZ, Cl, Cm, Cn = mapreduce(sum, vcat, surface_coefficients(system; axes = Wind()))
 
    return NF_COEFFS(CX, CY, CZ, Cl, Cm, Cn)
end

"""
    farfield_forces(system :: VortexLatticeSystem)

Compute the **farfield** forces in **wind axes** for all components of the `VortexLatticeSystem`.
"""
function farfield_forces(system :: VortexLatticeSystem)
    hs = system.vortices 
    Γs = system.circulations
    α  = system.freestream.alpha
    β  = system.freestream.beta
    V  = system.reference.speed
    ρ  = system.reference.density
    
    @views NamedTuple(key => farfield_forces(Γs[key], hs[key], V, α, β, ρ) for key in keys(hs))
end

"""
    farfield_coefficients(system :: VortexLatticeSystem)

Compute the **total farfield** force coefficients for all components of the `VortexLatticeSystem`. These are in **wind axes** by definition.
"""
farfield_coefficients(system :: VortexLatticeSystem) = let ; 
    map(farfield_forces(system)) do ff
        FF_COEFFS(force_coefficient(ff, dynamic_pressure(system.reference.density, system.reference.speed), system.reference.area))
    end
end

"""
    farfield(system :: VortexLatticeSystem)

Compute the **total farfield** force coefficients for all components of the `VortexLatticeSystem`. These are in **wind axes** by definition.
"""
function farfield(system :: VortexLatticeSystem)
    q = dynamic_pressure(system.reference.density, system.reference.speed)
    coeffs = force_coefficient(sum(farfield_forces(system)), q, system.reference.area)

    return FF_COEFFS(coeffs...)
end

function center_of_pressure(system :: VortexLatticeSystem)
    x_AC = system.reference.location[1]
    c_ref = system.reference.chord
    nf = nearfield(system)

    x_CP = x_AC - c_ref * nf.Cm / nf.CZ

    return x_CP
end

residual!(R, Γ, system :: VortexLatticeSystem) = solve_nonlinear!(R, system.vortices, Γ, -velocity(system.freestream), system.freestream.omega)