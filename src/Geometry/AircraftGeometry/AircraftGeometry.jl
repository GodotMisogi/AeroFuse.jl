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
using Setfield
using MacroTools

# Math tools
import ..MathTools: uniform_spacing, linear_spacing, sine_spacing, cosine_spacing, cosine_interp, splitat, adj3, slope, columns, forward_sum, forward_division, forward_difference, weighted_vector, vectarray, extend_yz

# Panel geometry
import ..PanelGeometry: Panel2D, Panel3D, panel_area, normal_vector, transform, make_panels, wetted_area

import ..AeroMDAO: properties

## Types
#==========================================================================================#

abstract type AbstractAircraft end

abstract type AbstractWing <: AbstractAircraft end

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

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
y_mac(y, b, λ) = y + b / 2 * (1 + 2λ) / (3(1 + λ))
quarter_chord(chord) = 0.25 * chord

include("Wings/halfwing.jl")
include("Wings/wing.jl")
include("Wings/mesh_tools.jl")
include("Wings/mesh_wing.jl")
include("Wings/controls.jl")
include("Wings/sections.jl")

"""
    aspect_ratio(wing :: AbstractWing)

Compute the aspect ratio of an `AbstractWing`.
"""
aspect_ratio(wing) = aspect_ratio(span(wing), projected_area(wing))

"""
    properties(wing :: AbstractWing)

Compute the generic properties of interest (span, area, etc.) of an `AbstractWing`.
"""
properties(wing :: AbstractWing) = [ aspect_ratio(wing), span(wing), projected_area(wing), mean_aerodynamic_chord(wing), mean_aerodynamic_center(wing) ]

function Base.show(io :: IO, wing :: AbstractWing)
    println(io, supertype(typeof(wing)),  " with ", length(spans(wing)), " spanwise section(s).")
    println(io, "Aspect Ratio: ", aspect_ratio(wing))
    println(io, "Span (m): ", span(wing))
    println(io, "Projected Area (m): ", projected_area(wing))
    println(io, "Mean Aerodynamic Chord (m): ", mean_aerodynamic_chord(wing))
    println(io, "Mean Aerodynamic Center (m): ", mean_aerodynamic_center(wing))

    nothing
end

end