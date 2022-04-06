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
y_mac(y, b, λ) = y + b / 2 * (1 + 2λ) / 3(1 + λ)
quarter_chord(chord) = 0.25 * chord

include("Wings/halfwing.jl")
include("Wings/wing.jl")
include("Wings/mesh_tools.jl")
include("Wings/mesh_wing.jl")
include("Wings/controls.jl")
include("Wings/sections.jl")

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
    taper_ratio(wing :: AbstractWing)

Compute the taper ratio of an `AbstractWing`, defined as the tip chord length divided by the root chord length.
"""
taper_ratio(wing :: AbstractWing) = taper_ratio(wing)

"""
    properties(wing :: AbstractWing)

Compute the generic properties of interest (span, area, etc.) of an `AbstractWing`.
"""
properties(wing :: AbstractWing) = [ aspect_ratio(wing), span(wing), projected_area(wing), mean_aerodynamic_chord(wing), mean_aerodynamic_center(wing) ]

"""
    mesh_chords(wing :: AbstractWing, n_s :: Vector{Integer}, n_c :: Integer; flip = false)

Mesh the span and chord distributions of an `AbstractWing` with ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions.
"""
mesh_chords(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_chords(wing, span_num, chord_num; spacings = spacings)

"""
    mesh_wing(wing :: AbstractWing, n_s :: Vector{Integer}, n_c :: Integer; flip = false)

Mesh the span and airfoil coordinate distributions of an `AbstractWing` with ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions.
"""
mesh_wing(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_wing(wing, span_num, chord_num; spacings = spacings)

"""
    mesh_cambers(wing :: AbstractWing, n_s :: Integer, n_c :: Integer; spacings = symmetric_spacing(wing))

Mesh the camber distribution of a `Wing` into panels of ``n_s`` spanwise divisions per section and ``n_c`` chordwise divisions with an `AbstractSpacing`` distribution.
"""
mesh_cambers(wing :: AbstractWing, span_num, chord_num; spacings = symmetric_spacing(wing)) = mesh_cambers(wing, span_num, chord_num; spacings = spacings)

"""
    coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the planform coordinates of a `HalfWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels..
"""
coordinates(wing :: AbstractWing, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine()) = chop_wing(coordinates(wing), span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

"""
    chord_coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the chord coordinates of an `AbstractWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
chord_coordinates(wing :: AbstractWing, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine()) = chord_coordinates(wing, span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

"""
    camber_coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the camber coordinates of an `AbstractWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
camber_coordinates(wing :: AbstractWing, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine()) = camber_coordinates(wing, span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

"""
    surface_coordinates(wing :: AbstractWing, n_s :: Integer, n_c :: Integer, flip = false)

Compute the surface coordinates of an `AbstractWing` given numbers of spanwise ``n_s`` and chordwise ``n_c`` panels, with an option to flip the signs of the ``y``-coordinates.
"""
surface_coordinates(wing :: AbstractWing, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine()) = surface_coordinates(wing, span_num, chord_num; span_spacing = span_spacing, chord_spacing = chord_spacing)

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