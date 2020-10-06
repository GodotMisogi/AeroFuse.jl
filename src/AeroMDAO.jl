module AeroMDAO

include("FoilParametrization.jl")
using LinearAlgebra
import Base: *, +
using Base.Iterators
using Statistics
using .FoilParametrization: Foil

#--------------------HACKS----------------------#

"""
"Lenses" to access subfields on lists of objects.
"""
|>(obj, fields :: Array{Symbol}) = foldl(getproperty, fields, init = obj)
|>(list_objs :: Array{T}, fields :: Array{Symbol}) where {T <: Any} = list_objs .|> [fields]

# Difference operators on lists?
fwdsum(xs) = xs[2:end] .+ xs[1:end-1]
fwddiff(xs) = xs[2:end] .- xs[1:end-1]
cendiff(xs) = xs[3:end] .- 2 * xs[2:end-1] .+ xs[1:end-2] 

#----------------VECTOR SPACES?---------------#

# Convert struct entries to lists
structtolist(x) = [ getproperty(x, name) for name ∈ (fieldnames ∘ typeof)(x) ]

abstract type Point end

*(scale :: Real, point :: Point) = typeof(point)(scale * structtolist(point)...)
+(p1 :: Point, p2 :: Point) = typeof(p1)(structtolist(p1) .+ structtolist(p2)...)

struct Point2D <: Point
    x :: Real; y :: Real;
end

struct Point3D <: Point
    x :: Real; y :: Real; z :: Real;
end

#----------------PANEL METHODS------------------#

abstract type Panels end

struct Panel <: Panels
    loc :: Array{Point}
end

"""
Evaluates the midpoint of the coordinates of a panel.
"""
# midpoint(panel :: Panel) = 


abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points, a chord length, and twist angle in degrees.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

#-----------------WING---------------------#

"""
Definition for a wing consisting of airfoils, span lengths, dihedrals, and sweep angles.
"""
struct HalfWing <: Aircraft
    foils :: Array{Foil} # Airfoil profiles
    chords :: Array{<: Real} # Chord lengths (m)
    spans :: Array{<: Real}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{<: Real} # Dihedral angles (deg)
    sweeps :: Array{<: Real} # Leading-edge sweep angles (deg)
    twists :: Array{<: Real} # Twist angles (deg)
    HalfWing(foils, chords, spans, dihedrals, sweeps, twists) = new(foils, chords, spans, 
    deg2rad.(dihedrals), deg2rad.(sweeps), deg2rad.(twists)) # Convert to radians
end

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(chord, taper_ratio) = (2/3) * chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
quarter_chord(chord) = 0.25 * chord

"""
Computes the span of a half-wing.
"""
span(wing :: HalfWing) = sum(wing.spans .* cos.(wing.dihedrals) .* cos.(wing.twists))

"""
Computes the projected area of a half-wing.
"""
function projected_area(wing :: HalfWing)
    mean_chords =  fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_spans = fwdsum(wing.spans) / 2 # Mean span lengths of sections
    sum(mean_chords .* mean_spans)
end

"""
Computes the mean aerodynamic chord of a half-wing.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2
    quarter_chords = (quarter_chord ∘ abs).(fwddiff(wing.chords))
    taper_ratios = map(/, wing.chords[2:end], wing.chords)
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord(wing.chords)
    sum(macs .* areas) / sum(areas)
end

struct Wing <: Aircraft
    left :: HalfWing
    Right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right
)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Wing) = aspect_ratio(span(wing), projected_area(wing))

# Kleisli-ish composition to transport global coordinates?
WingPos = Pair(Wing, Point3D)

end