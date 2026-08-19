## References
#==========================================================================================#

abstract type AbstractReferences end

"""
    References(V, ρ, μ, S, b, c, r)
    References(; 
        speed, density, viscosity,
        sound_speed, area, span, 
        chord, location
    )

Define reference values with speed ``V``, density ``ρ``, dynamic viscosity ``μ``, area ``S``, span ``b``, chord ``c``, location ``r`` for a vortex lattice analysis. A constructor with named arguments is provided for convenience:

# Arguments
- `speed       :: Real         = 1.`: Speed (m/s)
- `density     :: Real         = 1.225`: Density (m)
- `viscosity   :: Real         = 1.5e-5`: Dynamic viscosity (kg/(m ⋅ s))
- `sound_speed :: Real         = 330.`: Speed of sound (m/s)
- `area        :: Real         = 1.`: Area (m²)
- `span        :: Real         = 1.`: Span length (m)
- `chord       :: Real         = 1.`: Chord length (m)
- `location    :: Vector{Real} = [0,0,0]`: Position (m)
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
    return References{T}(V, ρ, μ, a, S, b, c, r) 
end

# References(; V, rho, mu, a, S, b, c, r) = References(V, rho, mu, a, S, b, c, r)

References(; speed = 1., density = 1.225, viscosity = 1.5e-5, sound_speed = 330., area = 1., span = 1., chord = 1., location = zeros(3)) = References(speed, density, viscosity, sound_speed, area, span, chord, location)

Base.broadcastable(refs :: References) = Ref(refs)

dynamic_pressure(refs :: References) = 1/2 * refs.density * refs.speed^2

kinematic_viscosity(refs :: References) = refs.viscosity / refs.density
mach_number(refs :: References) = refs.speed / refs.sound_speed
reynolds_number(refs :: References) = refs.density * refs.speed * refs.chord / refs.viscosity

force_coefficient(force, refs :: References) = force_coefficient(force, dynamic_pressure(refs), refs.area)
moment_coefficient(moment, refs :: References) = moment_coefficient(moment, dynamic_pressure(refs), refs.area, refs.span, refs.chord)

rate_coefficient(fs :: Freestream, refs :: References) = rate_coefficient(fs.omega, refs.speed, refs.span, refs.chord)

## System
#==========================================================================================#

abstract type AbstractPotentialFlowSystem end

"""
    VortexLatticeSystem

A system consisting of the relevant variables in a vortex lattice analysis for post-processing.

# Arguments
The accessible fields are:
- `vortices`:  The array of vortices, presently of `AbstractVortex` types.
- `circulations`: The circulation strengths of the vortices obtained by solving the linear system.
- `influence_matrix`: The influence matrix of the linear system.
- `boundary_vector`: The boundary condition corresponding to the right-hand-side of the linear system.
- `freestream :: Freestream`: The freestream conditions.
- `reference :: References`: The reference values.
"""
struct VortexLatticeSystem{
    M <: DenseArray{<: AbstractVortex},
    N,
    R,
    S,
    P,
    Q} <: AbstractPotentialFlowSystem
    vortices          :: M
    circulations      :: N 
    influence_matrix  :: R
    boundary_vector   :: S
    freestream        :: P
    reference         :: Q
    compressible      :: Bool
end

"""
    VortexLatticeSystem(
        aircraft, 
        fs :: Freestream, 
        refs :: References, 
        compressible = false, 
    )

Construct a `VortexLatticeSystem` for analyzing inviscid aerodynamics of an aircraft (must be a `ComponentArray` of `Horseshoe`s or `VortexRing`s) with `Freestream` conditions and `References` for non-dimensionalization. Options are provided for compressibility corrections via the Prandtl-Glauert transformation (`false` by default) and axis system for computing velocities and forces (`Geometry` by default).
"""
# Whether an assembled aircraft carries a fuselage line singularity block (`:fuse`).
has_fuselage(ac) = :fuse in propertynames(ac)

function VortexLatticeSystem(aircraft, fs :: Freestream, refs :: References, compressible = false, warn = true)

    M = mach_number(refs) # For Mach number bound checks

    # Compressible mode
    if compressible
        @assert M < 1. "Only compressible subsonic flow conditions (M < 1) are valid!"
        if M > 0.7 && warn
            # @warn "Results in transonic to sonic flow conditions (0.7 < M < 1) are most likely incorrect!" 
        end

        # (Prandtl-Glauert ∘ Wind axis) transformation
        β_pg = √(1 - M^2)
        ac = @. prandtl_glauert_scale_coordinates(geometry_to_wind_axes(aircraft, fs), β_pg)
    else # Incompressible mode
        if warn 
            # if M > 0.3 @warn "Compressible regime (M > 0.3) but compressibility correction is off, be wary of the analysis!" end
        end

        β_pg = 1
        ac = @. geometry_to_wind_axes(aircraft, fs)
    end

    # Quasi-steady freestream velocity
    U = geometry_to_wind_axes(-velocity(fs), fs.alpha, fs.beta)
    Ω = geometry_to_wind_axes(fs.omega, fs) / refs.speed

    # Solve system. If a fuselage line (`:fuse` block) is present, inject its prescribed
    # source (thickness) field into every collocation point's boundary condition, and solve
    # the fuselage doublet strengths together with the lifting-surface circulations.
    if has_fuselage(ac)
        Ups = map(el -> source_line_velocity(control_point(el), ac.fuse, U), ac)
        Γs, AIC, boco = solve_linear_fuselage(ac, U, Ups, Ω)
    else
        Γs, AIC, boco = solve_linear(ac, U, Ω)
    end

    return VortexLatticeSystem(aircraft, refs.speed * Γs / β_pg^2, AIC, boco, fs, refs, compressible)
end

# Miscellaneous
rate_coefficient(system :: VortexLatticeSystem) = rate_coefficient(system.freestream, system.reference)

## THINK ABOUT USING ONLY WIND AXES FOR PG-TRANSFORMATION AND MAPPING BACK

## Velocities
"""
    surface_velocities(
        system :: VortexLatticeSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the induced velocities for all components of the `VortexLatticeSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `VortexLatticeSystem` by default.
"""
function surface_velocities(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    vels = surface_velocities(system.vortices, system.vortices, system.circulations, system.reference.speed * -velocity(system.freestream), system.freestream.omega)
    return _vector_to_axes.(vels, Ref(axes), α, β)
end

## Forces
"""
    surface_forces(
        system :: VortexLatticeSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the forces for all components of the `VortexLatticeSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `VortexLatticeSystem` by default.
"""
function surface_forces(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    forces = surface_forces(system.vortices, system.circulations, system.reference.speed * -velocity(system.freestream), system.freestream.omega, system.reference.density)
    return _vector_to_axes.(forces, Ref(axes), α, β)
end

## Moments
"""
    surface_moments(
        system :: VortexLatticeSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the moments for all components of the `VortexLatticeSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `VortexLatticeSystem` by default.
"""
function surface_moments(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    geo_forces = surface_forces(system.vortices, system.circulations, system.reference.speed * -velocity(system.freestream), system.freestream.omega, system.reference.density)
    moments = surface_moments(system.vortices, geo_forces, system.reference.location)
    return _moment_to_axes.(moments, Ref(axes), α, β)
end

## Dynamics
"""
    surface_dynamics(
        system :: VortexLatticeSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the forces and moments for all components of the `VortexLatticeSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `VortexLatticeSystem` by default.
"""
function surface_dynamics(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    # Compute surface forces and moments in geometry axes
    surf_forces = surface_forces(system.vortices, system.circulations, system.reference.speed * -velocity(system.freestream), system.freestream.omega, system.reference.density)
    # The fuselage line carries no Kutta–Joukowsky force; substitute its slender-body
    # (Munk) sectional forces so the fuselage's own lift and destabilizing pitching moment
    # enter the nearfield, centre of pressure and stability derivatives.
    if has_fuselage(system.vortices)
        surf_forces.fuse .= fuselage_munk_forces(system)
    end
    surf_moments = surface_moments(system.vortices, surf_forces, system.reference.location)
    # Transform to target axes
    return _vector_to_axes.(surf_forces, Ref(axes), α, β), _moment_to_axes.(surf_moments, Ref(axes), α, β)
end

"""
    fuselage_munk_forces(system :: VortexLatticeSystem)

Slender-body (Munk) sectional force on each `FuselageLine` segment, in geometry axes. From
apparent-mass theory the sectional normal force is ``N'(x) = -\\tfrac{1}{2}\\rho U_\\infty\\,
d\\kappa/dx``, where ``\\kappa`` is the segment's cross-flow doublet strength (its stored,
dimensional circulation). Integrated over a closed body this gives zero net lift but a
destabilizing pitching moment.
"""
function fuselage_munk_forces(system :: VortexLatticeSystem)
    fuse = system.vortices.fuse
    κ    = system.circulations.fuse
    ρ    = system.reference.density
    U    = system.reference.speed
    n    = length(fuse)
    xs   = [ el.rc[1] for el in fuse ] # Axial stations (geometry axes)

    return map(1:n) do i
        # dκ/dx by finite difference (one-sided at the ends)
        dκdx = i == 1 ? (κ[2] - κ[1]) / (xs[2] - xs[1]) :
               i == n ? (κ[n] - κ[n-1]) / (xs[n] - xs[n-1]) :
                        (κ[i+1] - κ[i-1]) / (xs[i+1] - xs[i-1])
        Nʹ = -ρ * U / 2 * dκdx * segment_length(fuse[i])
        Nʹ * fuse[i].normal # Vertical (cross-flow) force in geometry axes
    end
end

"""
    surface_coefficients(
        system :: VortexLatticeSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the force and moment coefficients of the surfaces over all components in a given `VortexLatticeSystem`, in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `VortexLatticeSystem` by default.
"""
function surface_coefficients(system :: VortexLatticeSystem; axes :: AbstractAxisSystem = Geometry()) 
    # Compute surface forces in whichever axes
    forces, moments = surface_dynamics(system; axes)
    refs = system.reference

    # Compute coefficients
    CFs = @. force_coefficient(forces, refs)
    CMs = @. moment_coefficient(moments, dynamic_pressure(refs), refs.area, refs.span, refs.chord)

    return CFs, CMs
end

const NF_COEFFS = @SLArray (6) (:CX,:CY,:CZ,:Cl,:Cm,:Cn)
const FF_COEFFS = @SLArray (3) (:CDi,:CY,:CL)

"""
    nearfield_coefficients(system :: VortexLatticeSystem)

Compute the nearfield force and moment coefficients for all components of the `VortexLatticeSystem`. These are in **wind axes** by default.
"""
@views function nearfield_coefficients(system :: VortexLatticeSystem)
    # Compute surface force and moment coefficients in wind axes
    CFs, CMs = surface_coefficients(system; axes = Wind())
 
    # Construct NamedTuple with ComponentArray keys for each component
    return NamedTuple(key => NF_COEFFS(sum(CFs[key])..., sum(CMs[key])...) for key in keys(CFs))
end

"""
    nearfield(system :: VortexLatticeSystem)

Compute the **total** nearfield force and moment coefficients for all components of the `VortexLatticeSystem`. These are in **wind axes** by default.
"""
function nearfield(system :: VortexLatticeSystem)
    CFs, CMs = surface_coefficients(system; axes = Wind())
    return NF_COEFFS(vcat(sum(CFs), sum(CMs)))
end


"""
    farfield_forces(system :: VortexLatticeSystem)

Compute the **farfield** forces in **wind axes** for all components of the `VortexLatticeSystem`.
"""
@views function farfield_forces(system :: VortexLatticeSystem)
    hs = system.vortices 
    Γs = system.circulations
    α  = system.freestream.alpha
    β  = system.freestream.beta
    V  = system.reference.speed
    ρ  = system.reference.density
    
    # Construct NamedTuple with ComponentArray keys for each component. The fuselage line
    # (`:fuse`) has no trailing wake, so it contributes zero in the Trefftz plane — its
    # farfield force is handled by slender-body integration, not here. The key is kept (with a
    # zero force) so nearfield and farfield share the same component keys downstream.
    return NamedTuple(
        key => key == :fuse ? zero(SVector{3, typeof(V)}) : farfield_forces(Γs[key], hs[key], V, α, β, ρ)
        for key in keys(hs)
    )
end

"""
    farfield_coefficients(system :: VortexLatticeSystem)

Compute the **total farfield** force coefficients for all components of the `VortexLatticeSystem`. These are in **wind axes** by definition.
"""
farfield_coefficients(system :: VortexLatticeSystem) = map(farfield_forces(system)) do ff
        FF_COEFFS(force_coefficient(ff, system.reference))
    end

"""
    farfield(system :: VortexLatticeSystem)

Compute the **total farfield** force coefficients of the `VortexLatticeSystem`. These are in **wind axes** by definition.
"""
farfield(system :: VortexLatticeSystem) = FF_COEFFS(force_coefficient(sum(farfield_forces(system)), system.reference))

"""
    center_of_pressure(system :: VortexLatticeSystem)

Determine the center of pressure ``x_{cp}`` of the `VortexLatticeSystem`. 

This is computed based on the nearfield lift ``C_L`` and moment ``Cₘ`` coefficients, and the reference location ``xᵣ`` and chord length ``cᵣ`` from `References`: ``x_{cp} = xᵣ -cᵣ(Cₘ / C_L)``
"""
function center_of_pressure(system :: VortexLatticeSystem)
    x_ref = system.reference.location[1]
    c_ref = system.reference.chord
    nf = nearfield(system)

    x_CP = x_ref - c_ref * nf.Cm / nf.CZ

    return x_CP
end

# Consider adding spanwise loading later

# Residual equation for nonlinear analysis
residual!(R, Γ, system :: VortexLatticeSystem) = solve_nonlinear!(R, system.vortices, Γ, -velocity(system.freestream), system.freestream.omega)