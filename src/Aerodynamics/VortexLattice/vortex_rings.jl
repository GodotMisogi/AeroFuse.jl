"""
    VortexRing{T <: Real}(left_leg, bound_leg, back_leg, right_leg)

A vortex ring type consisting of four finite vortex line filaments connected in a loop.
"""
struct VortexRing{T <: Real} <: AbstractVortexArray
    left_leg  :: Line{T}
    bound_leg :: Line{T}
    back_leg  :: Line{T}
    right_leg :: Line{T}
end

"""
Helper function to compute the vortex ring given four points following Panel3D ordering.
"""
function vortex_ring(p1, p2, p3, p4)
    v1 = quarter_point(p1, p2)
    v4 = quarter_point(p4, p3)
    v2 = v1 + p2 - p1
    v3 = v4 + p3 - p4

    v1, v2, v3, v4
end

"""
Computes the vortex rings on a Panel3D.
"""
vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Constructor for vortex rings on a Panel3D using Lines. The following convention is adopted:

```
    p1 —bound_leg→ p4
    |               |
left_leg       right_leg
    ↓               ↓
    p2 —back_leg→ p3
```
"""
function VortexRing(panel :: Panel3D{T}) where T <: Real
    v1, v2, v3, v4 = vortex_ring(panel)

    bound_leg = Line(v1, v4)
    left_leg  = Line(v2, v1)
    right_leg = Line(v3, v2)
    back_leg  = Line(v4, v3)

    VortexRing{T}(left_leg, bound_leg, back_leg, right_leg)
end

Base.length(::VortexRing) = 1

"""
    sum_vortices(r, vortex_lines :: Array{Line}, Γ)

Sums the velocities evaluated at a point ``r`` of Lines with constant strength ``Γ``.
"""
sum_vortices(r, vortex_lines :: Array{<: Line}, Γ) = sum(line -> bound_velocity(r, line, Γ), vortex_lines)

"""
Computes the induced velocities at a point ``r`` of a VortexRing with constant strength ``Γ``.
"""
velocity(r, vortex_ring :: VortexRing, Γ) = sum_vortices(r, structtolist(vortex_ring), Γ)