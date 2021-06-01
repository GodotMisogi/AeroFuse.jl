module PanelGeometry

import Base: +, -, zero

using StaticArrays
using LinearAlgebra
using Rotations
using CoordinateTransformations

using ..AeroMDAO: Point2D, affine_2D, rotation, inverse_rotation, structtolist, sine_dist, cosine_dist, partition

export panel_dist, Panel, Panel2D, WakePanel2D, point1, point2, point3, point4, zero, collocation_point, panel_length, transform_panel, transform_panel_points, panel_angle, panel_tangent, panel_normal, panel_location, Panel3D, panel_area, panel_coords, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector

## Panel setup
#==========================================================================================#


abstract type Panel end

panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) - collocation_point(panel_1))

# split_panels(panels :: Vector{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

## 2D Panels
#==========================================================================================#

struct Panel2D{T <: Real} <: Panel
	p1 :: SVector{2,T}
	p2 :: SVector{2,T}
end

struct WakePanel2D{T <: Real} <: Panel
	p1 :: SVector{2,T}
	p2 :: SVector{2,T}
end

Panel2D(p1 :: FieldVector{2,T}, p2 :: FieldVector{2,T}) where T <: Real = Panel2D{T}(p1, p2)
WakePanel2D(p1 :: FieldVector{2,T}, p2 :: FieldVector{2,T}) where T <: Real = WakePanel2D{T}(p1, p2)

# Panel2D(p1 :: Union{Point2D{T}, MVector{2,T}, SVector{2,T}}, p2 :: Union{Point2D{T}, MVector{2,T}, SVector{2,T}}) where T <: Real = Panel2D{T}(Point2D(p1), Point2D(p2))

# Methods on panels
point1(p :: Panel) = p.p1
point2(p :: Panel) = p.p2
point3(p :: Panel) = p.p3
point4(p :: Panel) = p.p4

a :: Panel2D + b :: Panel2D = Panel2D(point1(a) + point1(b), point2(a) + point2(b))
a :: Panel2D - b :: Panel2D = Panel2D(point1(a) - point1(b), point2(a) - point2(b))

collocation_point(panel :: Panel, a = 0.5) = a * (point1(panel) + point2(panel))
panel_vector(panel :: Panel) = point2(panel) - point1(panel)
panel_length(panel :: Panel) = (norm ∘ panel_vector)(panel)

function transform_panel_points(panel_1 :: Panel, panel_2 :: Panel)
	x1, y1, x2, y2 = point1(panel_2), point2(panel_2)
	xs, ys = point1(panel_1)

	xp1, yp1 = affine_2D(x1, y1, xs, ys, panel_angle(panel_1)) 
	xp2, yp2 = affine_2D(x2, y2, xs, ys, panel_angle(panel_1)) 

	xp1, yp1, xp2, yp2
end

function transform_panel(panel :: Panel, point :: SVector{2,<: Real})
	xs, ys = point1(panel)
	affine_2D(first(point), last(point), xs, ys, panel_angle(panel)) 
end

transform_panel(panel_1 :: Panel, panel_2 :: Panel) = transform_panel(panel_1, collocation_point(panel_2))

panel_angle(panel :: Panel) = let (x, y) = panel_vector(panel); atan(y, x) end
panel_tangent(panel :: Panel) = rotation(1., 0., -panel_angle(panel))
panel_normal(panel :: Panel) = inverse_rotation(0., 1., panel_angle(panel))
panel_location(panel :: Panel) = let angle = panel_angle(panel); ifelse((π/2 <= angle <= π) || (-π <= angle <= -π/2), "lower", "upper") end

panel_points(panels :: Vector{<: Panel}) = [ point1.(panels); [(point2 ∘ last)(panels)] ]

reverse_panel(panel :: Panel) = Panel2D(panel.p2, panel.p1)

trailing_edge_panel(panels :: Vector{<: Panel}) = Panel2D((point2 ∘ last)(panels), (point1 ∘ first)(panels))

function wake_panel(panels, bound, α)
	firstx, firsty   = (point1 ∘ first)(panels)
	lastx, lasty	 = (point2 ∘ last)(panels)
	y_mid 			 = (firsty + lasty) / 2
	y_bound, x_bound = bound .* sincos(α) 
	WakePanel2D(SVector(lastx, y_mid), SVector(x_bound * lastx, y_bound * y_mid))
end

function wake_panels(panels, chord, length, num)
	firstx, firsty  = (point1 ∘ first)(panels)
	lastx, lasty	= (point2 ∘ last)(panels)
	y_mid 			= (firsty + lasty) / 2
	bounds 			= cosine_dist(chord + length / 2, length, num)
	@. WakePanel2D(SVector(bounds[1:end-1], y_mid), SVector(bounds[2:end], y_mid))
end

function panel_scalar(scalar_func, strength, panel :: Panel, x, y)
	# Transform point to local panel coordinates
    xp, yp = transform_panel(panel, SVector(x, y))
	scalar_func(strength, xp, yp, 0., panel_length(panel))
end

function panel_scalar(scalar_func, strength, panel_j :: Panel, panel_i :: Panel)
	x, y = collocation_point(panel_i)
	panel_scalar(scalar_func, strength, panel_j, x, y)
end

function panel_velocity(velocity_func, strength, panel :: Panel, x, y)
	# Transform point to local panel coordinates
    xp, yp = transform_panel(panel, SVector(x, y))

	# Compute velocity in local panel coordinates
	u, w = velocity_func(strength, xp, yp, 0., panel_length(panel))

	# Transform velocity to original coordinates
	inverse_rotation(u, w, panel_angle(panel))
end

function panel_velocity(velocity_func, strength, panel_j :: Panel, panel_i :: Panel, point = 0.5)
	x, y = collocation_point(panel_i, point)
	u, w = panel_velocity(velocity_func, strength, panel_j, x, y)
end

panel_velocity(f1, f2, str1, str2, panel :: Panel2D, x, y) = panel_velocity(f1, str1, panel, x, y) .+ panel_velocity(f2, str2, panel, x, y)

panel_velocity(f1, f2, str_j1, str_j2, panel_j :: Panel2D, panel_i :: Panel2D) = panel_velocity(f1, str_j1, panel_j, panel_i) .+ panel_velocity(f2, str_j2, panel_j, panel_i)

get_surface_values(panels, vals, surf = "upper", points = false) = partition(x -> (panel_location ∘ first)(x) == surf, (collect ∘ zip)(panels, vals), x -> ((first ∘ ifelse(points, point1, collocation_point) ∘ first)(x), last(x)))

## 3D Panels
#==========================================================================================#

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
struct Panel3D{T <: Real} <: Panel
	p1 :: SVector{3,T}
	p2 :: SVector{3,T}
	p3 :: SVector{3,T}
	p4 :: SVector{3,T}
end

Panel3D(p1 :: FieldVector{3,T}, p2 :: FieldVector{3,T}, p3 :: FieldVector{3,T}, p4 :: FieldVector{3,T}) where T <: Real = Panel3D{T}(p1, p2, p3, p4)

"""
	panel_area(panel :: Panel3D)

Compute the area of a Panel3D.
"""
panel_area(panel :: Panel3D) = norm(point2(panel) - point1(panel)) * norm(point3(panel) - point2(panel))

"""
	panel_coords(panel :: Panel3D)

Compute the coordinates of a Panel3D.
"""
panel_coords(panel :: Panel3D) = structtolist(panel)

"""
	transform(panel :: Panel3D, rotation, translation)

Perform an affine transformation on the coordinates of a Panel3D given a rotation matrix and translation vector.
"""
function transform(panel :: Panel3D, rotation, translation)
	p1, p2, p3, p4 = (Translation(translation) ∘ LinearMap(rotation)).(panel_coords(panel))	
	Panel3D(p1, p2, p3, p4)
end

"""
	midpoint(panel :: Panel3D)

Compute the midpoint of a Panel3D.
"""
midpoint(panel :: Panel3D) = (point1(panel) + point2(panel) + point3(panel) + point4(panel)) / 4

"""
	panel_normal(panel :: Panel3D)

Compute the normal vector of a Panel3D.
"""
panel_normal(panel :: Panel3D) = let p31 = point3(panel) - point1(panel), p42 = point4(panel) - point2(panel); p31 × p42 end

"""
	transform_normal(panel :: Panel3D, h_l, g_l)

Compute the normal ``n_l``, the normal ``n_0`` of a `Panel3D` perturbed by the control gain ``\\delta_l`` about the hinge axis ``h_l``.
"""
transform_normal(panel :: Panel3D, h_l, g_l) = g_l * cross(h_l, panel_normal(panel))
	
end
