module AircraftGeometry

using Base.Math
using Base.Iterators
using DelimitedFiles
using StaticArrays
using CoordinateTransformations
using Rotations
using LinearAlgebra

using ..AeroMDAO: uniform_spacing, linear_spacing, sine_spacing, cosine_spacing, cosine_interp, splitat, adj3, slope, columns, fwdsum, fwddiv, fwddiff, weighted_vector, vectarray, Point2D, Panel2D, Panel3D, extend_yz, transform, panel_area, panel_normal, wetted_area

abstract type Aircraft end

export Aircraft

## Foil geometry
#==========================================================================================#

include("foil.jl")

# export Foil, kulfan_CST, naca4, camber_CST, paneller, read_foil, split_foil, foil_camthick, camthick_foil, camber_thickness, cosine_foil, camthick_to_CST, coords_to_CST, max_thickness_to_chord_ratio_location

## Fuselage geometry

include("fuselage.jl")

# export Fuselage, projected_area, length, cosine_distribution

## Wing geometry
#==========================================================================================#

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
y_mac(y, b, λ) = y + b / 2 * (1 + 2λ) / 3(1 + λ)
quarter_chord(chord) = 0.25 * chord

include("halfwing.jl")
include("wing.jl")
include("mesh_tools.jl")
include("mesh_wing.jl")

aspect_ratio(wing) = aspect_ratio(span(wing), projected_area(wing))

info(wing :: Union{Wing, HalfWing}) = [ span(wing), projected_area(wing), mean_aerodynamic_chord(wing), aspect_ratio(wing) ]

# export HalfWing, HalfWingSection, Wing, WingSection, mean_aerodynamic_chord, span, aspect_ratio, projected_area, taper_ratio, info, max_tbyc_sweeps, leading_edge, trailing_edge, wing_bounds, chop_leading_edge, chop_trailing_edge, chop_wing, paneller, panel_wing, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing, mean_aerodynamic_center, wetted_area, number_of_spanwise_panels, spanwise_spacing, coordinates

end