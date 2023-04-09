# Velocity kernels
#==========================================================================================#

bound_leg_velocity(a, b, Γ)    = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))
trailing_leg_velocity(r, Γ, u) = Γ/4π * normalize(r) × normalize(u) / (norm(r) - dot(r, u))

trailing_legs_velocities(a, b, Γ, u) = trailing_leg_velocity(a, Γ, u) - trailing_leg_velocity(b, Γ, u)
total_horseshoe_velocity(a, b, Γ, u) = bound_leg_velocity(a, b, Γ) + trailing_legs_velocities(a, b, Γ, u)

# Finite-core velocity kernels
function bound_leg_velocity(a, b, Γ, ε) 
    na, nb, σ = norm(a), norm(b), dot(a, b)
    term_1 = (na^2 - σ) / √(na^2 + ε^2) + (nb^2 - σ) / √(nb^2 + ε^2)
    term_2 = a × b / (na^2 * nb^2 - σ^2 + ε^2 * (na^2 + nb^2 - 2 * na * nb))
    
    Γ/4π * term_1 * term_2
end

trailing_leg_velocity(r, Γ, u, ε) = Γ/4π * normalize(r) × u / (norm(r) - dot(r, u) + ε^2 / (norm(r) + dot(r, u)))
trailing_legs_velocities(a, b, Γ, u, ε) = trailing_leg_velocity(a, Γ, u, ε) - trailing_leg_velocity(b, Γ, u, ε)
total_horseshoe_velocity(a, b, Γ, u, ε) = bound_leg_velocity(a, b, Γ, ε) + trailing_legs_velocities(a, b, Γ, u, ε)

## Arrays of vortex lines
#==========================================================================================#

abstract type AbstractVortex end

## Horseshoe type
#==========================================================================================#

"""
    Horseshoe(r1, r2, rc, normal, chord)

Define a horseshoe vortex with a start and endpoints ``r₁, r₂`` for the bound leg, a collocation point ``r``, a normal vector ``n̂``, and a finite core size.

The finite core setup is not implemented for now.
"""
struct Horseshoe{T <: Real} <: AbstractVortex
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
    rc :: SVector{3,T}
    normal :: SVector{3,T}
    core :: T
end

Base.length(::Horseshoe) = 1

r1(hs :: Horseshoe) = hs.r1
r2(hs :: Horseshoe) = hs.r2

function Horseshoe(r1, r2, rc, n, c)
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(n), eltype(c))
    Horseshoe{T}(r1, r2, rc, n, c)
end

control_point(hs :: AbstractVortex) = hs.rc
normal_vector(hs :: AbstractVortex) = hs.normal

r1(r, hs :: Horseshoe) = r - hs.r1
r2(r, hs :: Horseshoe) = r - hs.r2

"""
    transform(hs :: Horseshoe, T :: LinearMap)

Generate a new `Horseshoe` with the points and normal vectors transformed by the `LinearMap` ``T``.
"""
transform(hs :: Horseshoe, T :: LinearMap) = setproperties(hs,
    r1 = T(hs.r1),
    r2 = T(hs.r2),
    rc = T(hs.rc),
    normal = T(hs.normal),
)

transform(hs :: Horseshoe; rotation = I(3), translation = zeros(3)) = transform(hs, Translation(translation) ∘ LinearMap(rotation))

"""
    bound_leg_center(hs :: Horseshoe)

Compute the midpoint of the bound leg of a `Horseshoe`.
"""
bound_leg_center(hs :: Horseshoe) = (hs.r1 + hs.r2) / 2

"""
    bound_leg_vector(hs :: Horseshoe)

Compute the direction vector of the bound leg of a `Horseshoe`.
"""
bound_leg_vector(hs :: Horseshoe) = hs.r2 - hs.r1

"""
    velocity(r, hs :: Horseshoe, Γ, u_hat = [1.,0.,0.])

Compute the induced velocity at a point ``r`` of a given `Horseshoe` with a bound leg of constant strength ``Γ`` and semi-infinite trailing legs pointing in a given direction ``û``, by default `û = x̂`.
"""
velocity(r, hs :: Horseshoe, Γ, V_hat = SVector{3,promote_type(eltype(r),eltype(Γ),eltype(hs.core))}(1,0,0)) = total_horseshoe_velocity(r - hs.r1, r - hs.r2, Γ, V_hat, hs.core)

"""
    bound_velocity(r, hs :: Horseshoe, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the bound leg with constant strength ``Γ`` of a given `Horseshoe`.
"""
bound_velocity(r, hs :: Horseshoe, Γ) = bound_leg_velocity(r - hs.r1, r - hs.r2, Γ, hs.core)

"""
    trailing_velocity(r, hs :: Horseshoes, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the semi-infinite trailing legs with constant strength ``Γ`` of a given `Horseshoe` `hs`.
"""
trailing_velocity(r, hs :: Horseshoe, Γ, V) = trailing_legs_velocities(r - hs.r1, r - hs.r2, Γ, V, hs.core)

## Vortex ring type
#==========================================================================================#

"""
    VortexRing(r1, r2, r3, r4, r_c, n̂, ε)

A vortex ring consisting of four points ``r_i, i = 1,…,4``, a collocation point ``r_c``, a normal vector ``n̂``, and a core size ``ε``. The following convention is adopted:

```
    r1 —front leg→ r4
    |               |
left leg       right leg
    ↓               ↓
    r2 —back leg-→ r3
```
"""
struct VortexRing{T <: Real} <: AbstractVortex
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
    r3 :: SVector{3,T}
    r4 :: SVector{3,T}
    rc :: SVector{3,T}
    normal :: SVector{3,T}
    trailing :: Bool
    core :: T
end

function VortexRing(r1, r2, r3, r4, rc, n, trailing, c)
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(n), eltype(c))
    VortexRing{T}(r1, r2, r3, r4, rc, n, trailing, c)
end

control_point(ring :: VortexRing) = ring.rc
normal_vector(ring :: VortexRing) = ring.normal

Base.length(::VortexRing) = 1

"""
    velocity(r, ring :: VortexRing, Γ)

Computes the velocity at a point ``r`` induced by a `VortexRing` with constant strength ``Γ``.
"""
function velocity(r, ring :: VortexRing, Γ, V_hat = SVector(1,0,0))
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

trailing_velocity(r, ring :: VortexRing, Γ, V) = trailing_legs_velocities(r - ring.r1, r - ring.r4, Γ, V, ring.core)

"""
    transform(ring :: VortexRing, T :: LinearMap)

Generate a new `VortexRing` with the points and normal vectors transformed by the `LinearMap` ``T``.
"""
transform(ring :: VortexRing, T :: LinearMap) = setproperties(ring,
    r1 = T(ring.r1),
    r2 = T(ring.r2),
    r3 = T(ring.r3),
    r4 = T(ring.r4),
    rc = T(ring.rc),
    normal = T(ring.normal),
    trailing = ring.trailing,
    core = ring.core,
)

bound_leg_center(ring :: VortexRing) = (ring.r1 + ring.r4) / 2
bound_leg_vector(ring :: VortexRing) = ring.r4 - ring.r1

r1(ring :: VortexRing) = ring.r1
r2(ring :: VortexRing) = ring.r4