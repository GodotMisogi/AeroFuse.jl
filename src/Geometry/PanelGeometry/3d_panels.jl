## 3D Panels
#==========================================================================================#

abstract type AbstractPanel3D <: AbstractPanel end

"""
    Panel3D(p1, p2, p3, p4)

A composite type consisting of 4 Cartesian coordinates `p1, p2, p3, p4` representing corners of a panel in 3 dimensions. The following commutative diagram depicts the order:

```
z → y
↓
x
        p1 —→— p4
        |       |
        ↓       ↓
        |       |
        p2 —→— p3
```
"""
struct Panel3D{T <: Real} <: AbstractPanel3D
    p1 :: SVector{3,T}
    p2 :: SVector{3,T}
    p3 :: SVector{3,T}
    p4 :: SVector{3,T}
end

Panel3D(p1 :: FieldVector{3,T}, p2 :: FieldVector{3,T}, p3 :: FieldVector{3,T}, p4 :: FieldVector{3,T}) where T <: Real = Panel3D{T}(p1, p2, p3, p4)

Panel3D((p1, p2, p3, p4)) = Panel3D(p1, p2, p3, p4)

Base.length(:: Panel3D) = 1

average_chord(panel :: Panel3D) = (p2(panel) - p1(panel) + p3(panel) - p4(panel)) / 2
average_width(panel :: Panel3D) = (p4(panel) - p1(panel) + p3(panel) - p2(panel)) / 2

"""
    panel_area(panel :: Panel3D)

Compute the area of a Panel3D.
"""
panel_area(panel :: Panel3D) = (norm ∘ cross)(average_chord(panel), average_width(panel))

wetted_area(panels :: Matrix{<: Panel3D}) = sum(panel -> panel_area(panel), panels)

"""
    panel_coords(panel :: Panel3D)

Compute the coordinates of a Panel3D.
"""
panel_coords(panel :: Panel3D) = [ p1(panel), p2(panel), p3(panel), p4(panel) ]

"""
    make_panels(xyzs)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
make_panels(xyzs) = @views Panel3D.(xyzs[1:end-1,1:end-1], xyzs[2:end,1:end-1], xyzs[2:end,2:end], xyzs[1:end-1,2:end])

"""
    transform(panel :: Panel3D, rotation, translation)

Perform an affine transformation on the coordinates of a Panel3D given a rotation matrix and translation vector.
"""
transform(panel :: Panel3D, rotation, translation) = Panel3D((Translation(translation) ∘ LinearMap(rotation)).(panel_coords(panel)))

"""
    midpoint(panel :: Panel3D)

Compute the midpoint of a Panel3D.
"""
midpoint(panel :: Panel3D) = (p1(panel) + p2(panel) + p3(panel) + p4(panel)) / 4

"""
    panel_normal(panel :: Panel3D)

Compute the normal vector of a Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = p3(panel) - p1(panel), p42 = p4(panel) - p2(panel); p31 × p42 end

"""
    transform_normal(panel :: Panel3D, h_l, g_l)

Compute the normal ``n_l``, the normal ``n_0`` of a `Panel3D` perturbed by the control gain ``\\delta_l`` about the hinge axis ``h_l``.
"""
transform_normal(panel :: Panel3D, h_l, g_l) = g_l * cross(h_l, panel_normal(panel))