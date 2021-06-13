module PanelGeometry

import Base: +, -, zero

using StaticArrays
using LinearAlgebra
using Rotations
using CoordinateTransformations

using ..AeroMDAO: Point2D, affine_2D, rotation, inverse_rotation, structtolist, sine_spacing, cosine_spacing, partition

export AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, AbstractPanel3D, Panel3D, p1, p2, p3, p4, zero, collocation_point, panel_length, transform_panel, transform_panel_points, panel_angle, panel_tangent, panel_normal, panel_location, panel_area, panel_coords, transform, midpoint, panel_points, wake_panel, wake_panels, reverse_panel, panel_velocity, panel_scalar, trailing_edge_panel, get_surface_values, panel_vector, panel_dist

## Panel setup
#==========================================================================================#


abstract type AbstractPanel end

include("2d_panels.jl")
include("3d_panels.jl")

# Methods on panels
p1(p :: AbstractPanel) = p.p1
p2(p :: AbstractPanel) = p.p2
p3(p :: AbstractPanel) = p.p3
p4(p :: AbstractPanel) = p.p4

transform_panel(panel_1 :: AbstractPanel, panel_2 :: AbstractPanel) = transform_panel(panel_1, collocation_point(panel_2))

panel_dist(panel_1 :: AbstractPanel, panel_2 :: AbstractPanel) = norm(collocation_point(panel_2) - collocation_point(panel_1))

# split_panels(panels :: Vector{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))
	
end
