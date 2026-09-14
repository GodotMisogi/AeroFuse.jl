# Velocity kernels
#==========================================================================================#

# For ModelingToolkit.jl/Symbolics.jl support with StaticArrays.jl norm method (see https://github.com/JuliaSymbolics/Symbolics.jl/issues/888)
# norm(v) = sqrt(sum(abs2, v))
# normalize(v) = v / norm(v)

bound_leg_velocity(a, b, Γ) = Γ / 4π * (1 / norm(a) + 1 / norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))
trailing_leg_velocity(r, Γ, u) = Γ / 4π * normalize(r) × normalize(u) / (norm(r) - dot(r, u))

trailing_legs_velocities(a, b, Γ, u) = trailing_leg_velocity(a, Γ, u) - trailing_leg_velocity(b, Γ, u)
total_horseshoe_velocity(a, b, Γ, u) = bound_leg_velocity(a, b, Γ) + trailing_legs_velocities(a, b, Γ, u)

# Finite-core velocity kernels
function bound_leg_velocity(a, b, Γ, ε)
    na, nb, σ = norm(a), norm(b), dot(a, b)
    term_1 = (na^2 - σ) / √(na^2 + ε^2) + (nb^2 - σ) / √(nb^2 + ε^2)
    term_2 = a × b / (na^2 * nb^2 - σ^2 + ε^2 * (na^2 + nb^2 - 2 * na * nb))

    Γ / 4π * term_1 * term_2
end

trailing_leg_velocity(r, Γ, u, ε) = Γ / 4π * normalize(r) × u / (norm(r) - dot(r, u) + ε^2 / (norm(r) + dot(r, u)))
trailing_legs_velocities(a, b, Γ, u, ε) = trailing_leg_velocity(a, Γ, u, ε) - trailing_leg_velocity(b, Γ, u, ε)
total_horseshoe_velocity(a, b, Γ, u, ε) = bound_leg_velocity(a, b, Γ, ε) + trailing_legs_velocities(a, b, Γ, u, ε)

## Arrays of vortex lines
#==========================================================================================#

abstract type AbstractVortex end

## Element models (solver-facing specifications)
#==========================================================================================#

"""
    AbstractElementModel

Supertype for the element-model specifications passed to [`elements`](@ref) to select the
singularity discretization for a component (e.g. [`Horseshoe`](@ref), [`VortexRing`](@ref),
[`SourcePanel`](@ref)). Each model is a lightweight, option-carrying tag; the populated
`AbstractVortex` elements it produces (`HorseshoeVortex`, `RingVortex`, `SourcePanel3D`) are
what the solver assembles into the influence system.
"""
abstract type AbstractElementModel end

## HorseshoeVortex type
#==========================================================================================#

"""
    Horseshoe(; core_size = 0.)

Element-model specification selecting horseshoe-vortex discretization. Pass it to
[`elements`](@ref) to build `HorseshoeVortex` elements from a mesh, e.g.
`elements(wing_mesh, Horseshoe())`. `core_size` sets the finite vortex-core radius.
"""
struct Horseshoe <: AbstractElementModel
    core_size::Float64
end

Horseshoe(; core_size = 0.) = Horseshoe(core_size)

"""
    HorseshoeVortex(r1, r2, rc, normal, chord)

Define a horseshoe vortex with a start and endpoints ``r₁, r₂`` for the bound leg, a collocation point ``r``, a normal vector ``n̂``, and a finite core size.

The finite core setup is not implemented for now.
"""
struct HorseshoeVortex{T} <: AbstractVortex
    r1::SVector{3,T}
    r2::SVector{3,T}
    rc::SVector{3,T}
    normal::SVector{3,T}
    core::T
end

Base.length(::HorseshoeVortex) = 1

r1(hs::HorseshoeVortex) = hs.r1
r2(hs::HorseshoeVortex) = hs.r2

function HorseshoeVortex(r1, r2, rc, n, c)
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(n), eltype(c))
    HorseshoeVortex{T}(r1, r2, rc, n, c)
end

control_point(hs::AbstractVortex) = hs.rc
normal_vector(hs::AbstractVortex) = hs.normal

"""
    has_wake(:: AbstractVortex)

Whether an element sheds a trailing wake and therefore contributes to the Trefftz-plane
farfield (induced-drag) integration. Lifting vortices shed a wake (`true`, the default);
non-lifting elements such as source panels and slender-body lines do not.
"""
has_wake(::AbstractVortex) = true

r1(r, hs::HorseshoeVortex) = r - hs.r1
r2(r, hs::HorseshoeVortex) = r - hs.r2

"""
    transform(hs :: HorseshoeVortex, T :: LinearMap)

Generate a new `HorseshoeVortex` with the points and normal vectors transformed by the `LinearMap` ``T``.
"""
transform(hs::HorseshoeVortex, T::LinearMap) = setproperties(hs,
    r1=T(hs.r1),
    r2=T(hs.r2),
    rc=T(hs.rc),
    normal=T(hs.normal),
)

transform(hs::HorseshoeVortex; rotation=I(3), translation=zeros(3)) = transform(hs, Translation(translation) ∘ LinearMap(rotation))

"""
    bound_leg_center(hs :: HorseshoeVortex)

Compute the midpoint of the bound leg of a `HorseshoeVortex`.
"""
bound_leg_center(hs::HorseshoeVortex) = (hs.r1 + hs.r2) / 2

"""
    bound_leg_vector(hs :: HorseshoeVortex)

Compute the direction vector of the bound leg of a `HorseshoeVortex`.
"""
bound_leg_vector(hs::HorseshoeVortex) = hs.r2 - hs.r1

"""
    velocity(r, hs :: HorseshoeVortex, Γ, u_hat = [1.,0.,0.])

Compute the induced velocity at a point ``r`` of a given `HorseshoeVortex` with a bound leg of constant strength ``Γ`` and semi-infinite trailing legs pointing in a given direction ``û``, by default `û = x̂`.
"""
velocity(r, hs::HorseshoeVortex, Γ, V_hat=SVector{3,promote_type(eltype(r), eltype(Γ), eltype(hs.core))}(1, 0, 0)) = total_horseshoe_velocity(r - hs.r1, r - hs.r2, Γ, V_hat, hs.core)

"""
    bound_velocity(r, hs :: HorseshoeVortex, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the bound leg with constant strength ``Γ`` of a given `HorseshoeVortex`.
"""
bound_velocity(r, hs::HorseshoeVortex, Γ) = bound_leg_velocity(r - hs.r1, r - hs.r2, Γ, hs.core)

"""
    trailing_velocity(r, hs :: HorseshoeVortexs, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the semi-infinite trailing legs with constant strength ``Γ`` of a given `HorseshoeVortex` `hs`.
"""
trailing_velocity(r, hs::HorseshoeVortex, Γ, V) = trailing_legs_velocities(r - hs.r1, r - hs.r2, Γ, V, hs.core)

## Vortex ring type
#==========================================================================================#

"""
    VortexRing(; core_size = 0., trailing = :auto)

Element-model specification selecting vortex-ring (lattice) discretization. Pass it to
[`elements`](@ref) to build `RingVortex` elements from a mesh, e.g.
`elements(wing_mesh, VortexRing())`. `core_size` sets the finite vortex-core radius;
`trailing` controls trailing-edge wake identification (`:auto` detects the trailing row of
the mesh, matching the historical behaviour).
"""
struct VortexRing <: AbstractElementModel
    core_size::Float64
    trailing::Symbol
end

VortexRing(; core_size = 0., trailing = :auto) = VortexRing(core_size, trailing)

"""
    RingVortex(r1, r2, r3, r4, r_c, n̂, ε)

A vortex ring consisting of four points ``r_i, i = 1,…,4``, a collocation point ``r_c``, a normal vector ``n̂``, and a core size ``ε``. The following convention is adopted:

```
    r1 —front leg→ r4
    |               |
left leg       right leg
    ↓               ↓
    r2 —back leg-→ r3
```
"""
struct RingVortex{T} <: AbstractVortex
    r1::SVector{3,T}
    r2::SVector{3,T}
    r3::SVector{3,T}
    r4::SVector{3,T}
    rc::SVector{3,T}
    normal::SVector{3,T}
    trailing::Bool
    core::T
end

function RingVortex(r1, r2, r3, r4, rc, n, trailing, c)
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(n), eltype(c))
    RingVortex{T}(r1, r2, r3, r4, rc, n, trailing, c)
end

control_point(ring::RingVortex) = ring.rc
normal_vector(ring::RingVortex) = ring.normal

Base.length(::RingVortex) = 1

"""
    velocity(r, ring :: RingVortex, Γ)

Computes the velocity at a point ``r`` induced by a `RingVortex` with constant strength ``Γ``.
"""
function velocity(r, ring::RingVortex, Γ, V_hat=SVector(1, 0, 0))
    # Compute vectors to evaluation point
    r1, r2, r3, r4 = r - ring.r1, r - ring.r2, r - ring.r3, r - ring.r4
    core = ring.core

    # Evaluate bound leg velocities and sum
    v1 = bound_leg_velocity(r1, r4, Γ, core)
    v2 = bound_leg_velocity(r4, r3, Γ, core)
    v3 = bound_leg_velocity(r3, r2, Γ, core)
    v4 = bound_leg_velocity(r2, r1, Γ, core)

    v = v1 + v2 + v3 + v4

    # Add horseshoe velocity contribution if trailing edge panel
    if ring.trailing
        v += total_horseshoe_velocity(r2, r3, Γ, V_hat, ring.core)
    end

    return v
end

trailing_velocity(r, ring::RingVortex, Γ, V) = trailing_legs_velocities(r - ring.r1, r - ring.r4, Γ, V, ring.core)

"""
    transform(ring :: RingVortex, T :: LinearMap)

Generate a new `RingVortex` with the points and normal vectors transformed by the `LinearMap` ``T``.
"""
transform(ring::RingVortex, T::LinearMap) = setproperties(ring,
    r1=T(ring.r1),
    r2=T(ring.r2),
    r3=T(ring.r3),
    r4=T(ring.r4),
    rc=T(ring.rc),
    normal=T(ring.normal),
    trailing=ring.trailing,
    core=ring.core,
)

bound_leg_center(ring::RingVortex) = (ring.r1 + ring.r4) / 2
bound_leg_vector(ring::RingVortex) = ring.r4 - ring.r1

r1(ring::RingVortex) = ring.r1
r2(ring::RingVortex) = ring.r4