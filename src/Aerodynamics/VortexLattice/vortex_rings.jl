"""
    VortexRing(r1, r2, r3, r4, r_c, n̂, ε)

A vortex ring consisting of four points ``r_i, i = 1,…,4``, a collocation point ``r_c``, a normal vector ``n̂``, and a core size ``ε``.
"""
struct VortexRing{T <: Real} <: AbstractVortex
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
    r3 :: SVector{3,T}
    r4 :: SVector{3,T}
    rc :: SVector{3,T}
    normal :: SVector{3,T}
    core :: T
end

control_point(ring :: VortexRing) = ring.rc
normal_vector(ring :: VortexRing) = ring.normal

function vortex_lines(p1, p2, p3, p4)
    v1 = quarter_point(p1, p2)
    v4 = quarter_point(p4, p3)
    v2 = v1 + p2 - p1
    v3 = v4 + p3 - p4

    v1, v2, v3, v4
end

"""
Constructor for vortex rings on a Panel3D using Lines. The following convention is adopted:

```
    p1 —front leg→ p4
    |               |
left leg       right leg
    ↓               ↓
    p2 —back leg-→ p3
```
"""
function VortexRing(panel :: Panel3D{T}, normal = normal_vector(panel), ε = 0.) where T <: Real
    r_c = control_point(panel)
    VortexRing{T}(panel.p1, panel.p2, panel.p3, panel.p4, r_c, normal, ε)
end

Base.length(::VortexRing) = 1

"""
Computes the induced velocities at a point ``r`` of a VortexRing with constant strength ``Γ``.
"""
function velocity(r, ring :: VortexRing, Γ) 
    v1 = bound_leg_velocity(r - ring.r1, r - ring.r2, Γ, ring.core)
    v2 = bound_leg_velocity(r - ring.r2, r - ring.r3, Γ, ring.core)
    v3 = bound_leg_velocity(r - ring.r3, r - ring.r4, Γ, ring.core)
    v4 = bound_leg_velocity(r - ring.r4, r - ring.r1, Γ, ring.core)

    return v1 + v2 + v3 + v4
end

# influence_matrix(rings) = [ influence_coefficient(ring_j, ring_i) for ring_j in rings, ring_i in rings ] 

# boundary_condition(rings, U, Ω) =  map(hs -> dot(U + Ω × control_point(hs), ring_normal(hs)), rings)