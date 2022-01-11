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

# Math tools
import ..MathTools: uniform_spacing, linear_spacing, sine_spacing, cosine_spacing, cosine_interp, splitat, adj3, slope, columns, fwdsum, fwddiv, fwddiff, weighted_vector, vectarray, extend_yz

# Panel geometry
import ..PanelGeometry: Panel2D, Panel3D, panel_area, panel_normal, transform, make_panels

import ..AeroMDAO: properties

## Types
#==========================================================================================#

abstract type AbstractAircraft end

abstract type AbstractWing <: AbstractAircraft end

## Foil geometry
#==========================================================================================#

include("foil.jl")

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
y_mac(y, b, λ) = y + b / 2 * (1 + 2λ) / 3(1 + λ)
quarter_chord(chord) = 0.25 * chord

include("halfwing.jl")
include("wing.jl")
include("mesh_tools.jl")
include("mesh_wing.jl")



"""
    span(wing :: AbstractWing)

Compute the planform span of an `AbstractWing`.
"""
span(wing :: AbstractWing) = span(wing)

"""
    projected_area(wing :: AbstractWing)

Compute the projected area of an `AbstractWing`` by summing the trapezoidal areas.
"""
projected_area(wing :: AbstractWing) = projected_area(wing)


"""
    mean_aerodynamic_chord(wing :: AbstractWing)

Compute the mean aerodynamic chord of an `AbstractWing`.
"""
mean_aerodynamic_chord(wing :: AbstractWing) = mean_aerodynamic_chord(wing)

"""
    mean_aerodynamic_chord(wing :: AbstractWing)

Compute the coordinates of the mean aerodynamic center of an `AbstractWing`.
"""
mean_aerodynamic_center(wing :: AbstractWing) = mean_aerodynamic_center(wing)

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
    println(io, "$(supertype(typeof(wing))) with $(length(spans(wing))) spanwise section(s).")
    AR, b, S, c, mac = properties(wing)
    println(io, "Aspect Ratio: $AR")
    println(io, "Span (m): $b")
    println(io, "Projected Area (m): $S")
    println(io, "Mean Aerodynamic Chord (m): $c")
    println(io, "Mean Aerodynamic Center (m): $mac")

    nothing
end

end