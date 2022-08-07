module PanelGeometry

## Package imports
#==========================================================================================#

import Base: +, -, zero

using StaticArrays
using LinearAlgebra
import LinearAlgebra: ×
using Rotations
using CoordinateTransformations
using SplitApplyCombine

import ..Laplace: AbstractFreestream, Freestream, velocity
import ..MathTools: affine_2D, rotation, inverse_rotation, structtolist, sine_spacing, cosine_spacing, partition, Point2D, Point3D

## Panel setup
#==========================================================================================#

abstract type AbstractPanel end

abstract type AbstractPanel3D <: AbstractPanel end

include("2d_panels.jl")
include("3d_panels.jl")

# Methods on panels
p1(p :: AbstractPanel) = p.p1
p2(p :: AbstractPanel) = p.p2
p3(p :: AbstractPanel3D) = p.p3
p4(p :: AbstractPanel3D) = p.p4

xs(p :: AbstractPanel3D) = SVector(p.p1.x, p.p2.x, p.p3.x, p.p4.x)
ys(p :: AbstractPanel3D) = SVector(p.p1.y, p.p2.y, p.p3.y, p.p4.y)
zs(p :: AbstractPanel3D) = SVector(p.p1.z, p.p2.z, p.p3.z, p.p4.z)

transform_panel(panel_1 :: AbstractPanel, panel_2 :: AbstractPanel) = transform_panel(panel_1, collocation_point(panel_2))

distance(panel_1 :: AbstractPanel, panel_2 :: AbstractPanel) = norm(collocation_point(panel_2) - collocation_point(panel_1))

# split_panels(panels :: Vector{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

end
