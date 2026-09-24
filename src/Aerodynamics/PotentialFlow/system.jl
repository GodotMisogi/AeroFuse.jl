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

Define reference values with speed ``V``, density ``ρ``, dynamic viscosity ``μ``, area ``S``, span ``b``, chord ``c``, location ``r`` for a potential-flow analysis. A constructor with named arguments is provided for convenience:

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
    PotentialFlowSystem

A coupled velocity-based potential-flow solution for aerodynamic post-processing.

# Arguments
The accessible fields are:
- `elements`: Named component arrays of populated `AbstractPotentialFlowElement`s.
- `strengths`: Solved strengths with the same component layout as `elements`.
  Horseshoes and rings store circulation, source panels store source density,
  and fuselage lines store cross-flow doublet-line strength.
- `influence_matrix`: The influence matrix of the linear system.
- `boundary_vector`: The boundary condition corresponding to the right-hand-side of the linear system.
- `freestream :: Freestream`: The freestream conditions.
- `reference :: References`: The reference values.
- `slipstream`: The prescribed propeller slipstream(s) for blown-lift analyses (`nothing` if absent).
- `backend`: The compute backend used for the solve and post-processing (`nothing` for the serial host path).
  With a device backend, `influence_matrix` and `boundary_vector` remain resident on that device.
"""
struct PotentialFlowSystem{
    M <: DenseArray{<: AbstractPotentialFlowElement},
    N,
    R,
    S,
    P,
    Q,
    W,
    B} <: AbstractPotentialFlowSystem
    elements          :: M
    strengths         :: N
    influence_matrix  :: R
    boundary_vector   :: S
    freestream        :: P
    reference         :: Q
    compressible      :: Bool
    slipstream        :: W
    backend           :: B
end

# Property aliases preserve existing source-level access without changing field layout.
@inline function Base.getproperty(system::PotentialFlowSystem, name::Symbol)
    field = name === :vortices ? :elements :
            name === :circulations ? :strengths : name
    return getfield(system, field)
end

Base.propertynames(system::PotentialFlowSystem, private::Bool = false) =
    (fieldnames(typeof(system))..., :vortices, :circulations)

function _validate_fuselage_stations(block, key)
    length(block) >= 2 || throw(ArgumentError(
        "Fuselage component $key requires at least two axial stations."))
    xs = map(el -> control_point(el)[1], block)
    all(isfinite, xs) && all(>(0), diff(vec(xs))) || throw(ArgumentError(
        "Fuselage component $key requires finite, strictly increasing axial stations."))
    return nothing
end

function _validate_components(aircraft)
    aircraft isa ComponentArray || throw(ArgumentError(
        "Aircraft must be a ComponentArray of named potential-flow element arrays."))
    isempty(aircraft) && throw(ArgumentError("Aircraft must contain at least one element."))
    for key in keys(aircraft)
        block = aircraft[key]
        block isa AbstractArray || throw(ArgumentError(
            "Component $key must be an array of populated potential-flow elements."))
        isempty(block) && continue
        element = first(block)
        element isa AbstractPotentialFlowElement || throw(ArgumentError(
            "Component $key must contain populated potential-flow elements, not model tags."))
        all(el -> typeof(el) === typeof(element), block) || throw(ArgumentError(
            "Component $key must contain one concrete element type; split mixed types into separate components."))
        element isa FuselageLine && _validate_fuselage_stations(block, key)
    end
    return nothing
end

"""
    PotentialFlowSystem(
        aircraft, 
        fs :: Freestream, 
        refs :: References, 
        compressible = false, 
    )

Solve inviscid aerodynamics for named `ComponentArray` components of populated
potential-flow elements, such as `HorseshoeVortex`, `RingVortex`, `SourcePanel3D`,
and `FuselageLine`. Each component contains one concrete element type. Fuselage
lines require at least two strictly increasing axial stations.

`compressible` enables the subsonic Prandtl–Glauert transformation; `warn`
controls regime warnings. `slipstream` supplies prescribed propeller disks.
Post-processing axes are chosen on the evaluation methods, not the constructor.
The current solver applies the same speed/PG strength scaling to all families;
compressible non-vortex strength recovery has not been independently validated.

`backend` selects where the influence matrix is assembled and solved and where the ``O(N^2)``
induced-velocity post-processing runs: `nothing` (default) is the serial host path; a
`KernelAbstractions.Backend` (e.g. `CPU()` for multithreading, `MetalBackend()`, `CUDABackend()`)
requires loading `KernelAbstractions` (or a GPU package providing it). Backends without Float64
support compute in Float32. Device backends do not support automatic differentiation.
"""
function PotentialFlowSystem(aircraft, fs :: Freestream, refs :: References, compressible = false, warn = true; slipstream = nothing, backend = nothing)

    _validate_components(aircraft)

    M = mach_number(refs) # For Mach number bound checks

    # Normalize the slipstream input to a collection of disks (or `nothing`) so it can be
    # broadcast/summed uniformly downstream.
    slips = slipstream isa PropellerDisk ? [slipstream] : slipstream

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

    # Solve system. Prescribed onset-flow fields are injected into every collocation point's
    # boundary condition (per unit freestream speed, in the solved wind/PG axes): the fuselage
    # thickness source (from any `FuselageLine` elements) and/or the propeller slipstream(s). The
    # slipstream disks are transformed into the same axes as `ac` so they can be evaluated at its
    # collocation points. Any element with a non-standard boundary condition (e.g. a
    # `FuselageLine`'s slender-body cylinder relation) rewrites its own row inside `solve_linear`,
    # so a single generic solve handles every combination — the block names carry no behaviour.
    has_slip   = !isnothing(slips)
    slips_ac   = has_slip ? map(d -> prandtl_glauert_scale_coordinates(geometry_to_wind_axes(d, fs.alpha, fs.beta), β_pg), slips) : slips
    fuse_elems = filter(el -> el isa FuselageLine, ac)   # prescribed thickness-source line, if any
    Ups = map(ac) do el
        cp = control_point(el)
        v_fuse = isempty(fuse_elems) ? zero(cp) : source_line_velocity(cp, fuse_elems, U)
        v_slip = has_slip ? slipstream_velocity(cp, slips_ac, refs) / refs.speed : zero(cp)
        v_fuse + v_slip
    end
    Γs, AIC, boco = device_solve_linear(backend, ac, U, Ups, Ω)

    return PotentialFlowSystem(aircraft, refs.speed * Γs / β_pg^2, AIC, boco, fs, refs, compressible, slips, backend)
end

# Miscellaneous
rate_coefficient(system :: PotentialFlowSystem) = rate_coefficient(system.freestream, system.reference)

## THINK ABOUT USING ONLY WIND AXES FOR PG-TRANSFORMATION AND MAPPING BACK

## Slipstream (blown lift)
# Per-panel propeller-slipstream velocity at the bound-leg centres (geometry axes), added to
# the nearfield Kutta–Joukowsky velocity for the local dynamic-pressure boost. Zero everywhere
# when the system carries no slipstream.
_slipstream_velocity(r, ::Nothing, refs) = zero(r)
_slipstream_velocity(r, slips, refs) = slipstream_velocity(r, slips, refs)

bound_slipstream_velocities(system :: PotentialFlowSystem) = map(el -> _slipstream_velocity(bound_leg_center(el), system.slipstream, system.reference), system.elements)

## Backend dispatch for O(N²) post-processing
# `nothing` keeps the serial host path verbatim. Device backends evaluate the induced sums via
# `device_induced_sum` and apply the same freestream/rotation/slipstream terms on the host.

# Copy device results (possibly lower precision) into the component layout of `like`.
_with_layout(like, values) = copyto!(similar(like), values)

_nearfield_velocities(::Nothing, system, U, Ω, Vps) = surface_velocities(system.elements, system.elements, system.strengths, U, Ω, Vps)

function _nearfield_velocities(backend, system, U, Ω, Vps)
    pts = map(bound_leg_center, system.elements)
    ind = _with_layout(pts, device_induced_sum(backend, trailing_velocity, pts, system.elements, system.strengths, -normalize(U)))
    return map((v, r, Vp) -> v - (U + Ω × r) + Vp, ind, pts, Vps)
end

_nearfield_forces(::Nothing, system, U, Ω, ρ, Vps) = surface_forces(system.elements, system.strengths, U, Ω, ρ, Vps)

function _nearfield_forces(backend, system, U, Ω, ρ, Vps)
    vels = _nearfield_velocities(backend, system, U, Ω, Vps)
    return map((h, Γ, V) -> kutta_joukowsky(ρ, V, bound_leg_vector(h), Γ), system.elements, system.strengths, vels)
end

## Velocities
"""
    surface_velocities(
        system :: PotentialFlowSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Evaluate freestream, rotation, trailing-induced flow, and propeller slipstream
at each element's bound-leg centre (control point for non-vortex adapters).
Bound-vortex induction and prescribed fuselage thickness sources are excluded.
Use `body_surface_velocities` for source-panel pressure evaluation.

The reference axis system is set to the geometry axes defined in the construction of the `PotentialFlowSystem` by default.
"""
function surface_velocities(system :: PotentialFlowSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    Vps = bound_slipstream_velocities(system)
    vels = _nearfield_velocities(system.backend, system, system.reference.speed * -velocity(system.freestream), system.freestream.omega, Vps)
    return _vector_to_axes.(vels, Ref(axes), α, β)
end

## Forces
"""
    surface_forces(
        system :: PotentialFlowSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the forces for all components of the `PotentialFlowSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `PotentialFlowSystem` by default.
"""
function surface_forces(system :: PotentialFlowSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    forces = _geometry_surface_forces(system)
    return _vector_to_axes.(forces, Ref(axes), α, β)
end

## Moments
"""
    surface_moments(
        system :: PotentialFlowSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the moments for all components of the `PotentialFlowSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `PotentialFlowSystem` by default.
"""
function surface_moments(system :: PotentialFlowSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    geo_forces = _geometry_surface_forces(system)
    moments = surface_moments(system.elements, geo_forces, system.reference.location)
    return _moment_to_axes.(moments, Ref(axes), α, β)
end

## Dynamics
"""
    surface_dynamics(
        system :: PotentialFlowSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the forces and moments for all components of the `PotentialFlowSystem` in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `PotentialFlowSystem` by default.
"""
function surface_dynamics(system :: PotentialFlowSystem; axes :: AbstractAxisSystem = Geometry())
    α, β = system.freestream.alpha, system.freestream.beta
    surf_forces = _geometry_surface_forces(system)
    surf_moments = surface_moments(system.elements, surf_forces, system.reference.location)
    # Transform to target axes
    return _vector_to_axes.(surf_forces, Ref(axes), α, β), _moment_to_axes.(surf_moments, Ref(axes), α, β)
end

function _geometry_surface_forces(system::PotentialFlowSystem)
    # Compute surface forces and moments in geometry axes
    Vps = bound_slipstream_velocities(system)
    surf_forces = _nearfield_forces(system.backend, system, system.reference.speed * -velocity(system.freestream), system.freestream.omega, system.reference.density, Vps)
    # Elements that carry no Kutta–Joukowsky force (source panels, slender-body lines) substitute
    # their own nearfield force model so their lift, pressure drag and pitching moment enter the
    # nearfield, centre of pressure and stability derivatives. The substitution is dispatched on
    # the block's element type (via `first`), not its name, so component names stay cosmetic.
    for key in propertynames(system.elements)
        block = getproperty(system.elements, key)
        isempty(block) && continue
        forces = block_force_override(first(block), system, key)
        isnothing(forces) || (getproperty(surf_forces, key) .= forces)
    end
    return surf_forces
end

"""
    block_force_override(element, system :: PotentialFlowSystem, key)

Nearfield surface force for the component block named `key`, when its elements do not obey the
Kutta–Joukowsky force model assumed by the generic nearfield loop. Dispatched on the block's
element type (`element`, the first element of the block) rather than the block name, so
component names carry no behaviour. Returns `nothing` for standard lifting vortices (the
Kutta–Joukowsky force already computed for the block stands), or the substituted per-element
forces (geometry axes) otherwise.
"""
block_force_override(::AbstractPotentialFlowElement, system :: PotentialFlowSystem, key) = nothing
block_force_override(::FuselageLine, system :: PotentialFlowSystem, key)   = fuselage_munk_forces(system, key)
block_force_override(::SourcePanel3D, system :: PotentialFlowSystem, key)  = body_forces(system, key)

"""
    fuselage_munk_forces(system :: PotentialFlowSystem, key = :fuse)

Slender-body (Munk) sectional force on each `FuselageLine` segment in the block `key`, in
geometry axes. From apparent-mass theory the sectional normal force is
``N'(x) = -\\tfrac{1}{2}\\rho U_\\infty\\, d\\kappa/dx``, where ``\\kappa`` is the segment's
stored dimensional cross-flow doublet-line strength. Integrated over a closed
body this gives zero net lift but a destabilizing pitching moment.
"""
function fuselage_munk_forces(system :: PotentialFlowSystem, key = :fuse)
    fuse = getproperty(system.elements, key)
    _validate_fuselage_stations(fuse, key)
    κ    = getproperty(system.strengths, key)
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
    body_surface_velocities(system :: PotentialFlowSystem, key = :body)

Flow velocity from freestream, rotation, and solved element strengths at the control point of every
`SourcePanel3D` in the block `key`, in geometry axes. Prescribed fuselage thickness sources and
propeller slipstreams are excluded. On a converged Neumann body the normal
component is ≈ 0, so this is essentially the tangential surface velocity.
"""
function body_surface_velocities(system :: PotentialFlowSystem, key = :body)
    U = system.reference.speed * -velocity(system.freestream)
    Ω = system.freestream.omega
    return _body_velocities(system.backend, system, getproperty(system.elements, key), U, Ω)
end

_body_velocities(::Nothing, system, block, U, Ω) = map(el -> induced_velocity(control_point(el), system.elements, system.strengths, U, Ω), block)

function _body_velocities(backend, system, block, U, Ω)
    pts = map(control_point, block)
    ind = _with_layout(pts, device_induced_sum(backend, velocity, pts, system.elements, system.strengths, -normalize(U)))
    return map((v, r) -> v - (U + Ω × r), ind, pts)
end

"""
    body_pressure_coefficients(system :: PotentialFlowSystem, key = :body)

Incompressible pressure coefficient ``C_p = 1 - (V_t/V_\\infty)^2`` at each body source panel
in the block `key`, from the tangential surface velocity ``V_t`` (see
[`body_surface_velocities`]).
"""
function body_pressure_coefficients(system :: PotentialFlowSystem, key = :body)
    V    = system.reference.speed
    vels = body_surface_velocities(system, key)
    return map(getproperty(system.elements, key), vels) do el, vel
        Vt² = dot(vel, vel) - dot(vel, normal_vector(el))^2
        1 - Vt² / V^2
    end
end

"""
    body_forces(system :: PotentialFlowSystem, key = :body)

Pressure force ``F = -C_p\\, q_\\infty A\\, n̂`` on each body source panel in the block `key`, in
geometry axes, from the surface pressure coefficients (see [`body_pressure_coefficients`]).
Summed over the closed body these give its lift, pressure drag and Munk pitching moment.
"""
function body_forces(system :: PotentialFlowSystem, key = :body)
    q   = dynamic_pressure(system.reference)
    Cps = body_pressure_coefficients(system, key)
    return map((el, Cp) -> -Cp * q * el.area * normal_vector(el), getproperty(system.elements, key), Cps)
end

"""
    body_forces(panels, elements, strengths, U, Ω, V, q)

Pressure force ``F = -C_p\\, q\\, A\\, n̂`` (geometry axes) on each body source panel in `panels`,
evaluated directly from a raw induced field rather than a solved `PotentialFlowSystem`: the
total velocity at each panel control point is induced by the full element list `elements` with
strengths `strengths` in the freestream `U`/`Ω`, giving ``C_p = 1 - V_t^2/V^2`` from the tangential
speed. `V` is the freestream speed and `q` the dynamic pressure.

This is the low-level twin of [`body_forces`](@ref)`(system)`, for coupled solvers that carry
the deformed lifting vortices and rigid body panels in one shared influence system.
"""
function body_forces(panels, elements, strengths, U, Ω, V, q)
    map(panels) do el
        vel = induced_velocity(control_point(el), elements, strengths, U, Ω)
        Vt² = dot(vel, vel) - dot(vel, normal_vector(el))^2
        Cp  = 1 - Vt² / V^2
        -Cp * q * el.area * normal_vector(el)
    end
end

"""
    surface_coefficients(
        system :: PotentialFlowSystem; 
        axes :: AbstractAxisSystem = Geometry()
    )

Compute the force and moment coefficients of the surfaces over all components in a given `PotentialFlowSystem`, in a specified reference axis system as a named argument.

The reference axis system is set to the geometry axes defined in the construction of the `PotentialFlowSystem` by default.
"""
function surface_coefficients(system :: PotentialFlowSystem; axes :: AbstractAxisSystem = Geometry()) 
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
    nearfield_coefficients(system :: PotentialFlowSystem)

Compute the nearfield force and moment coefficients for all components of the `PotentialFlowSystem`. These are in **wind axes** by default.
"""
@views function nearfield_coefficients(system :: PotentialFlowSystem)
    # Compute surface force and moment coefficients in wind axes
    CFs, CMs = surface_coefficients(system; axes = Wind())
 
    # Construct NamedTuple with ComponentArray keys for each component
    return NamedTuple(key => NF_COEFFS(sum(CFs[key])..., sum(CMs[key])...) for key in keys(CFs))
end

"""
    nearfield(system :: PotentialFlowSystem)

Compute the **total** nearfield force and moment coefficients for all components of the `PotentialFlowSystem`. These are in **wind axes** by default.
"""
function nearfield(system :: PotentialFlowSystem)
    CFs, CMs = surface_coefficients(system; axes = Wind())
    return NF_COEFFS(vcat(sum(CFs), sum(CMs)))
end


"""
    farfield_forces(system :: PotentialFlowSystem)

Compute the **farfield** forces in **wind axes** for all components of the `PotentialFlowSystem`.
"""
@views function farfield_forces(system :: PotentialFlowSystem)
    hs = system.elements 
    Γs = system.strengths
    α  = system.freestream.alpha
    β  = system.freestream.beta
    V  = system.reference.speed
    ρ  = system.reference.density
    
    # Construct NamedTuple with ComponentArray keys for each component. Elements with no
    # trailing wake (source panels, slender-body lines) contribute zero in the Trefftz plane —
    # their farfield force is handled separately (or is identically zero). The key is kept (with
    # a zero force) so nearfield and farfield share the same component keys downstream. The wake
    # test is dispatched on the block's element type (via `first`), not its name.
    return NamedTuple(
        key => (isempty(hs[key]) || !has_wake(first(hs[key]))) ? zero(SVector{3, typeof(V)}) : farfield_forces(Γs[key], hs[key], V, α, β, ρ)
        for key in keys(hs)
    )
end

"""
    farfield_coefficients(system :: PotentialFlowSystem)

Compute the **total farfield** force coefficients for all components of the `PotentialFlowSystem`. These are in **wind axes** by definition.
"""
farfield_coefficients(system :: PotentialFlowSystem) = map(farfield_forces(system)) do ff
        FF_COEFFS(force_coefficient(ff, system.reference))
    end

"""
    farfield(system :: PotentialFlowSystem)

Compute the **total farfield** force coefficients of the `PotentialFlowSystem`. These are in **wind axes** by definition.
"""
farfield(system :: PotentialFlowSystem) = FF_COEFFS(force_coefficient(sum(farfield_forces(system)), system.reference))

"""
    center_of_pressure(system :: PotentialFlowSystem)

Determine the center of pressure ``x_{cp}`` of the `PotentialFlowSystem`. 

This is computed based on the nearfield lift ``C_L`` and moment ``Cₘ`` coefficients, and the reference location ``xᵣ`` and chord length ``cᵣ`` from `References`: ``x_{cp} = xᵣ -cᵣ(Cₘ / C_L)``
"""
function center_of_pressure(system :: PotentialFlowSystem)
    x_ref = system.reference.location[1]
    c_ref = system.reference.chord
    nf = nearfield(system)

    x_CP = x_ref - c_ref * nf.Cm / nf.CZ

    return x_CP
end

# Consider adding spanwise loading later

# Residual equation for nonlinear analysis
residual!(R, Γ, system :: PotentialFlowSystem) = solve_nonlinear!(R, system.elements, Γ, -velocity(system.freestream), system.freestream.omega)