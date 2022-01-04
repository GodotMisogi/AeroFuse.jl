module PanelGeometry

## Package imports
#==========================================================================================#

import Base: +, -, zero

using StaticArrays
using LinearAlgebra
using Rotations
using CoordinateTransformations
using SplitApplyCombine

import ..MathTools: affine_2D, rotation, inverse_rotation, structtolist, sine_spacing, cosine_spacing, partition

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
