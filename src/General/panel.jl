module PanelGeometry

import Base: +, -, zero

using StaticArrays
using LinearAlgebra

include("../General/math_tools.jl")
using .MathTools: span, structtolist

## Panel setup
#==========================================================================================#

"""
Placeholder. Panels should be an abstract type as they have some common methods, at least in 2D and 3D. 
"""
abstract type Panel end

# Crap for automatic differentiation
# zero(:: NTuple{2,<:Real}) = (0., 0.)
# +(::Union{Nothing, Panel2D}, ::Union{Nothing,Panel2D}) = nothing
# zero(:: Nothing) = nothing

## 2D Panels
#==========================================================================================#

struct Panel2D <: Panel
    p1 :: NTuple{2, Real}
    p2 :: NTuple{2, Real}
end

# Methods on panels
point1(p :: Panel) = p.p1
point2(p :: Panel) = p.p2
zero(:: Panel2D) = Panel2D((0., 0.), (0., 0.))


a :: Panel + b :: Panel = Panel2D(point1(a) + point1(b), point2(a) + point2(b))
a :: Panel - b :: Panel = Panel2D(point1(a) - point1(b), point2(a) - point2(b))

collocation_point(panel :: Panel2D) = (point1(panel) .+ point2(panel)) ./ 2
panel_length(panel :: Panel2D) = norm(point2(panel) .- point1(panel))

panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) .- collocation_point(panel_1))
split_panels(panels :: Array{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

paneller(coords :: Array{<: Real, 2}) = [ Panel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ eachrow([ coords[2:end,:] coords[1:end-1,:] ]) ][end:-1:1]


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
struct Panel3D <: Panel
    p1 :: SVector{3, Real}
    p2 :: SVector{3, Real}
    p3 :: SVector{3, Real}
    p4 :: SVector{3, Real}
end

"""
Computes the area of Panel3D.
"""
panel_area(panel :: Panel3D) = norm(panel.p2 .- panel.p1) * norm(panel.p3 .- panel.p2)

"""
Computes the coordinates of a Panel3D.
"""
panel_coords(panel :: Panel3D) = structtolist(panel)

"""
Performs an affine transformation on the coordinates of a Panel3D.
"""
transform(panel :: Panel3D; rotation = one(RotMatrix{3, Float64}), translation = SVector(0, 0, 0)) = Panel3D(( (rotation * coord) .+ translation for coord in panel_coords(panel) )...)

"""
Computes the midpoint of Panel3D.
"""
midpoint(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) ./ 4

"""
Computes the normal vector of Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = panel.p3 .- panel.p1, 
                                     p42 = panel.p4 .- panel.p2, 
                                     p31_x_p42 = cross(p31, p42);
                                     p31_x_p42 end
                                     
end
