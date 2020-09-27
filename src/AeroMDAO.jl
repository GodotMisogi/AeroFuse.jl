module AeroMDAO

using LinearAlgebra
using Base.Iterators
using Interpolations
using Statistics

# Point definitions and methods

abstract type Point end

struct Point3D{T <: Real} <: Point
    x :: T
    y :: T
    z :: T
end

struct Point2D{T <: Real} <: Point
    x :: T
    y :: T
end

# Panel definitions and methods 

abstract type Panels end

struct Panel <: Panels
    loc :: Array{Point}
end

"""
Evaluates the midpoint of the coordinates of a panel.
"""
midpoint(panel :: Panel) = mean.(map(a -> getproperty.(panel.loc, a), (fieldnames ∘ typeof ∘ first)(panel.loc)))

# Aircraft parameters

abstract type Aircraft end

"""
Definition for a wing section consisting of its position, chord length, twist, and an airfoil.
"""
mutable struct WingSection <: Aircraft    
    loc :: Point # Coordinates of the section's leading edge with respect to the reference frame of the wing.
    chord :: Float64  # Chord length of the section.
    twist :: Float64  # Twist of the section about the leading edge.
    airfoil :: Array{Point} # The airfoil profile as an array of coordinates.
end

"""
Definition of a wing consisting of its position and sections.
"""
struct Wing <: Aircraft
    loc :: Point # Origin of the reference frame of the wing.
    sections :: Array{WingSection} # Collection of wing sections.
    # function Wing(ref)
    #     loc = ref
    #     # sections = [] ???
    #     return new(loc, sections)
end

"""
Computes the projected area of a wing given its sections.
"""
function projected_area(wing :: Wing)
    secs = wing.sections
    chords = getproperty.(secs, :chord)
    spans = getproperty.(getproperty.(secs, :loc), :y) 
    mean_chords = map(mean, chords, chords[:end-1]) # Mean chord lengths of sections.
    mean_spans = map(mean, spans, )
    return sum(mean_chords * mean_spans)
end

end