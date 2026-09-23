module AircraftGeometry

## Package imports
#==========================================================================================#

using Base.Math
using Base.Iterators
using DelimitedFiles
using StaticArrays
using CoordinateTransformations
using Rotations
using LinearAlgebra
using SplitApplyCombine
using Interpolations
using Accessors
using MacroTools

# Math tools
import ..MathTools: uniform_spacing, linear_spacing, sine_spacing, cosine_spacing, cosine_interp, splitat, adj3, slope, columns, forward_sum, forward_division, forward_difference, weighted_vector, vectarray, extend_yz

# Panel geometry
import ..PanelGeometry: Panel2D, Panel3D, panel_area, normal_vector, transform, make_panels, wetted_area

## Types
#==========================================================================================#

abstract type AbstractAircraft end

abstract type AbstractWing <: AbstractAircraft end
abstract type AbstractFuselage <: AbstractAircraft end

## Foil geometry
#==========================================================================================#

include("Foils/foil.jl")
include("Foils/class_shape_transformation.jl")
include("Foils/naca_airfoils.jl")

## Fuselage geometry
#==========================================================================================#

include("fuselage.jl")

## Wing geometry
#==========================================================================================#

# include("Wings/controls.jl")
include("Wings/halfwing.jl")
include("Wings/mesh_tools.jl")
include("Wings/mesh_wing.jl")
include("Wings/sections.jl")

end