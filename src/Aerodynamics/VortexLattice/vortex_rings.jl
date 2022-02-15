"""
    VortexRing(r1, r2, r3, r4, r_c, n̂, ε)

A vortex ring consisting of four points ``r_i, i = 1,…,4``, a collocation point ``r_c``, a normal vector ``n̂``, and a core size ``ε``.
"""
struct VortexRing{T <: Real} <: AbstractVortex
    r1                :: SVector{3,T}
    r2                :: SVector{3,T}
    r3                :: SVector{3,T}
    r4                :: SVector{3,T}
    collocation_point :: SVector{3,T}
    normal            :: SVector{3,T}
    core              :: T
end

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
function VortexRing(panel :: Panel3D{T}, normal = panel_normal(panel)) where T <: Real
    r_c = collocation_point(panel)
    ε   = 0. # (norm ∘ average_chord)(panel)
    VortexRing{T}(panel.p1, panel.p2, panel.p3, panel.p4, r_c, normal, ε)
end

Base.length(::VortexRing) = 1

"""
Computes the induced velocities at a point ``r`` of a VortexRing with constant strength ``Γ``.
"""
velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)