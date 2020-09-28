module AeroMDAO

using LinearAlgebra
import Base: *, +
using Base.Iterators
using Interpolations
using Statistics
using DelimitedFiles

# "Lenses" to access subfields on lists of objects
|>(obj, fields :: Array{Symbol}) = foldl(getproperty, fields, init = obj)
|>(list_objs :: Array{T}, fields :: Array{Symbol}) where {T <: Any} = list_objs .|> [fields]

# Convert tuple entries to lists
tupletolist(x) = [ getproperty(x, name) for name in (collect ∘ fieldnames)(typeof(x)) ]

#------------------POINT---------------#

abstract type Point end


struct Point2D{T <: Real} <: Point
    x :: T; y :: T;
end

*(scale :: Real, point :: Point2D) = Point2D(scale .* tupletolist(point)...)

struct Point3D{T <: Real} <: Point
    x :: T; y :: T; z :: T;
end

*(scale :: Real, point :: Point3D) = Point3D(scale .* tupletolist(point)...)


#----------------PANEL------------------#

abstract type Panels end

struct Panel <: Panels
    loc :: Array{Point}
end

"""
Evaluates the midpoint of the coordinates of a panel.
"""
midpoint(panel :: Panel) = mean.(map(a -> getproperty.(panel.loc, a), (fieldnames ∘ typeof ∘ first)(panel.loc)))


abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points, a chord length, and twist angle in degrees.
"""
struct Foil <: Aircraft
    foil :: Array{Point} # The foil profile as an array of coordinates, must be in Selig format.
    chord :: Float64  # Chord length of the section.
    twist :: Float64  # Twist of the profile about the leading edge (in degrees).
    Foil(foil, chord, twist = 0) = new(chord * foil, chord, deg2rad(twist))
end

"""
Reads a '.dat' file consisting of 2D coordinates, for an airfoil.
"""
function read_foil(path :: String)
    f = readdlm(path, skipstart = 1)
    Point2D{Float64}.(f[:,1], f[:,2])
end

"""
Provides the projections to the x-axis for a circle of given diameter and center.
"""
function cosine_dist(x_center :: Real, diameter :: Real, n :: Integer = 40)

    x_center .+ (diameter / 2.) .* cos.(range(-π, stop = 0, length = n))
end

"""
Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
"""
function cosine_foil(airfoil :: Foil, n :: Integer = 40)

    af = tupletolist.(airfoil.foil)
    #     first.(af), last.(af) end

    # d = maximum(x) - minimum(x)
    # x_center = (maximum(x) + minimum(x)) / 2

    # x_circ = cosine_dist(x_center, d, n)
    # itp_circ = LinearInterpolation(xs, ys)
    # y_circ = itp_circ(x_circ)

    # Panel(Point2D.(x_circ, y_circ))
end

#-----------------WING---------------------#

"""
Definition for a wing section consisting of root and tip airfoils, with the span length between them.
"""
struct WingSection <: Aircraft
    root :: Foil     # Root airfoil
    tip :: Foil      # Tip airfoil
    span :: Float64  # Leading edge to leading edge
end

fwdsum(xs, scale = 1) = scale * (xs[1:end-1] .+ xs[2:end])
aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(chord, taper_ratio) = (2/3) * chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)

"""
Definition of a wing consisting of its position and sections.
"""
struct Wing <: Aircraft
    sections :: Array{WingSection} # Collection of wing sections.
end

WingPos = Pair(Wing, Point3D)

"""
Computes the projected area of a wing given its sections.
"""
function projected_area(wing :: Wing)
    secs = wing.sections
    chords = secs |> [:chord]
    spans = secs |> [:loc, :y] 
    mean_chords = fwdsum(chords, 0.5) # Mean chord lengths of sections.
    mean_spans = fwdsum(spans, 0.5) # Mean span lengths of sections
    sum(mean_chords .* mean_spans)
end

function mean_aerodynamic_chord(wing :: Wing)
    secs = wing.sections
    chords = secs |> [:chord]
    mean_chords = fwdsum(chords, 0.5)
    quarter_chords = 0.25 * (chords[1:end-1] .- chords[2:end])
    spans = secs |> [:loc, :y] 
    area = sum(mean_chords .* spans)
end

end