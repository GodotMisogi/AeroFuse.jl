module AeroMDAO

include("MathTools.jl")
include("PanelMethods.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum
using .PanelMethods: Panel3D
using StaticArrays
using Rotations

abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

#-----------------WING---------------------#

"""
Definition for a half-wing consisting of airfoils, span lengths, dihedrals, and sweep angles.
"""
struct HalfWing <: Aircraft
    foils :: Array{Foil} # Airfoil profiles
    chords :: Array{<: Real} # Chord lengths (m)
    spans :: Array{<: Real}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{<: Real} # Dihedral angles (deg)
    sweeps :: Array{<: Real} # Leading-edge sweep angles (deg)
    twists :: Array{<: Real} # Twist angles (deg)
    HalfWing(foils, chords, spans, dihedrals, sweeps, twists) = new(foils, chords, spans, deg2rad.(dihedrals), deg2rad.(sweeps), -deg2rad.(twists)) # Convert to radians
end

aspect_ratio(span, area) = span^2 / area
mean_geometric_chord(span, area) = area / span
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord
area(span, chord) = span * chord
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
quarter_chord(chord) = 0.25 * chord

"""
Computes the planform span of a half-wing.
"""
span(wing :: HalfWing) = sum(wing.spans .* cos.(wing.dihedrals) .* cos.(fwdsum(wing.twists) / 2))

"""
Computes the projected area of a half-wing.
"""
function projected_area(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2 # Mean chord lengths of sections.
    mean_twists = fwdsum(wing.twists) / 2 # Mean twist angles of sections
    sum(@. wing.spans * cos(wing.dihedrals) * cos(wing.sweeps) * mean_chords * cos(mean_twists))
end

"""
Computes the mean aerodynamic chord of a half-wing.
"""
function mean_aerodynamic_chord(wing :: HalfWing)
    mean_chords = fwdsum(wing.chords) / 2
    taper_ratios = wing.chords[2:end] ./ wing.chords[1:end-1]
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord.(wing.chords[1:end-1], taper_ratios)
    sum(macs .* areas) / sum(areas)
end

function wing_coords(wing :: HalfWing)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps
    root_chord, tip_chords = peel(wing.chords)
    root_twist, tip_twists = peel(wing.twists)

    sweeped_spans, dihedraled_spans, cum_spans = cumsum(spans .* sin.(sweeps)), cumsum(spans .* sin.(dihedrals)), cumsum(spans)
    twisted_chords = tip_chords .* sin.(tip_twists)

    leading_xyz = [ sweeped_spans cum_spans dihedraled_spans ]
    trailing_xyz = [ (tip_chords .+ sweeped_spans) (cum_spans) (dihedraled_spans .+ twisted_chords) ]

    leading = [ 0 0 0; leading_xyz ]
    trailing = [ root_chord 0 root_chord * sin(root_twist); trailing_xyz ]

    SVector.(leading[:,1], leading[:,2], leading[:,3]),
    SVector.(trailing[:,1], trailing[:,2], trailing[:,3])
end

function wing_sections(wing :: HalfWing)
    lead, trail = wing_coords(wing)
    
    [ [le; te] for (le, te) in zip(lead, trail) ] 
end

function make_panels(wing :: HalfWing)
    lead, trail = wing_coords(wing);
    xs = [ lead[1:end-1] lead[2:end] trail[2:end] trail[1:end-1] ]

    [ Panel3D(pt...) for pt in eachrow(xs) ] 
end

# For horseshoe collocations
horseshoe_points(wing :: HalfWing) = [ 1/4 * fwdsum(wing.chords) / 2 cumsum(wing.spans / 2) ]
horseshoe_collocation(wing :: HalfWing) = [ 3/4 * fwdsum(wing.chords) /2 cumsum(wing.spans / 2) ]

struct Wing <: Aircraft
    left :: HalfWing
    right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right
)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Wing) = aspect_ratio(span(wing), projected_area(wing))

function wing_coords(wing :: Wing)
    left_lead, left_trail = wing_coords(wing.left)
    right_lead, right_trail = wing_coords(wing.right)

    leading = [ left_lead; right_lead ]
    trailing = [ left_trail; right_trail ]

    yflip!(leading), yflip!(trailing)
end

function yflip!(xs)
    print(xs) 
    xs[:,2] .= -xs[:,2]

    xs
end

wing_sections(wing :: Wing) = yflip!(wing_sections(wing.left)), wing_sections(wing.right)


# Kleisli-ish composition to transport global coordinates?
# PlanePos = Pair(Aircraft :: Aircraft, Point3D)

end