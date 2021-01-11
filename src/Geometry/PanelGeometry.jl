module PanelGeometry

import Base: +, -, zero

using StaticArrays
using LinearAlgebra
using Rotations
using CoordinateTransformations
using Zygote

using ..AeroMDAO: Point2D, affine_2D, rotation, inverse_rotation, span

export panel_dist, split_panels, Panel2D, point1, point2, point3, point4, zero, collocation_point, panel_length, trans_panel, panel_angle, panel_tangent, panel_normal, panel_location, Panel3D, panel_area, panel_coords, transform, midpoint

## Panel setup
#==========================================================================================#

"""
Placeholder. Panels should be an abstract type as they have some common methods, at least in 2D and 3D. 
"""
abstract type Panel end

panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) - collocation_point(panel_1))

split_panels(panels :: AbstractVector{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

## 2D Panels
#==========================================================================================#

struct Panel2D{T} <: Panel where T <: Real
    p1 :: Union{Point2D{T}, SVector{2,T}, MVector{2,T}}
    p2 :: Union{Point2D{T}, SVector{2,T}, MVector{2,T}}
end

Panel2D(p1 :: Union{Point2D{T}, SVector{2,T}, MVector{2,T}}, p2 :: Union{Point2D{T}, SVector{2,T}, MVector{2,T}}) where T <: Real = Panel2D{T}(p1, p2)

# Methods on panels
point1(p :: Panel) = p.p1
point2(p :: Panel) = p.p2
point3(p :: Panel) = p.p3
point4(p :: Panel) = p.p4

a :: Panel2D + b :: Panel2D = Panel2D(point1(a) + point1(b), point2(a) + point2(b))
a :: Panel2D - b :: Panel2D = Panel2D(point1(a) - point1(b), point2(a) - point2(b))

zero(::Panel2D) = Panel2D(Point2D(0.,0.), Point2D(0.,0.))
@Zygote.adjoint point1(p::Panel2D) = p.p1, x̄ -> (Panel2D(x̄, Point2D(0., 0.)),)
@Zygote.adjoint point2(p::Panel2D) = p.p2, ȳ -> (Panel2D(Point2D(0., 0.), ȳ),)
@Zygote.adjoint Panel2D(a, b) = Panel2D(a, b), p̄ -> (p̄.p1, p̄.p2)

collocation_point(panel :: Panel2D) = (point1(panel) + point2(panel)) / 2
panel_length(panel :: Panel2D) = norm(point2(panel) - point1(panel))

function trans_panel(panel_1 :: Panel2D, panel_2 :: Panel2D)
    x, y = collocation_point(panel_2)
    xs, ys = point1(panel_1)
    affine_2D(x, y, xs, ys, panel_angle(panel_1)) 
end

function panel_angle(panel :: Panel2D)
    xs, ys = point1(panel) 
    xe, ye = point2(panel)
    
    atan(ye - ys, xe - xs) 
end

panel_tangent(panel :: Panel2D) = rotation(1., 0., -panel_angle(panel))
panel_normal(panel :: Panel2D) = inverse_rotation(0., 1., panel_angle(panel))
panel_location(panel :: Panel2D) = let angle = panel_angle(panel); (π/2 <= angle <= π) || (-π <= angle <= -π/2) ? "lower" : "upper" end

## 3D Panels
#==========================================================================================#

"""
A composite type consisting of 4 coordinates. The following ASCII art depicts the order:

z → y
↓
x
        p1 —→— p4
        |       |
        ↓       ↓
        |       |
        p2 —→— p3
"""
struct Panel3D{T <: Real} <: Panel
    p1 :: SVector{3,T}
    p2 :: SVector{3,T}
    p3 :: SVector{3,T}
    p4 :: SVector{3,T}
end

Panel3D(p1 :: SVector{3,T}, p2 :: SVector{3,T}, p3 :: SVector{3,T}, p4 :: SVector{3,T}) where T <: Real = Panel3D{T}(p1, p2, p3, p4)

"""
Computes the area of Panel3D.
"""
panel_area(panel :: Panel3D) = norm(panel.p2 - panel.p1) * norm(panel.p3 - panel.p2)

"""
Computes the coordinates of a Panel3D.
"""
panel_coords(panel :: Panel3D) = structtolist(panel)

"""
Performs an affine transformation on the coordinates of a Panel3D.
"""
transform(panel :: Panel3D, rotation, translation) = Panel3D( (Translation(translation) ∘ LinearMap(rotation)).(panel_coords(panel))...)

"""
Computes the midpoint of Panel3D.
"""
midpoint(panel :: Panel3D) = (panel.p1 + panel.p2 + panel.p3 + panel.p4) / 4

"""
Computes the normal vector of Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = panel.p3 - panel.p1, 
                                     p42 = panel.p4 - panel.p2;
                                     p31 × p42 end
                                     
end
