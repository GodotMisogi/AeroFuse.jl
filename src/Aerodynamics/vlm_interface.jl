function quarter_point(p1, p2) 
    μ = SVector(1/4, 0, 1/4)
    return @. (1 - μ) * p1 + μ * p2
end 

function three_quarter_point(p1, p2) 
    μ = SVector(3/4, 0, 3/4)
    return @. (1 - μ) * p1 + μ * p2
end 

control_point(p1, p2, p3, p4) = (three_quarter_point(p1, p2) + three_quarter_point(p4, p3)) / 2
bound_leg(p1, p2, p3, p4) = (quarter_point(p1, p2), quarter_point(p4, p3))

"""
    bound_leg(panel :: Panel3D)

Compute the bound leg for a `Panel3D`, for horseshoes/vortex rings.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
    control_point(panel :: Panel3D)

Compute the control point of a `Panel3D` for horseshoes/vortex rings, which is the 3-quarter point on each side in the ``x``-``z`` plane.
"""
control_point(panel :: Panel3D) = control_point(panel.p1, panel.p2, panel.p3, panel.p4)


"""
    Horseshoe(panel :: Panel3D, normal, drift = zeros(3))

Generate a `Horseshoe` corresponding to a `Panel3D`, an associated normal vector, and a "drift velocity".
"""
function Horseshoe(panel :: Panel3D, normal, drift = SVector(0., 0., 0.); core_size = 0.)
    r1, r2 = bound_leg(panel)
    rc = control_point(panel) + drift
    Horseshoe(r1, r2, rc, normal, core_size)
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
    # r1 = quarter_point(panel.p1, panel.p2)
    # r4 = quarter_point(panel.p4, panel.p3)
    # r2 = normalize(panel.p2 - panel.p1) * 0.25 + r1
    # r3 = normalize(panel.p3 - panel.p4) * 0.25 + r2
    rc = control_point(panel) # (r1 + r2 + r3 + r4) / 4
    VortexRing{T}(panel.p1, panel.p2, panel.p3, panel.p4, rc, normal, ε)
end