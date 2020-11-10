module AeroMDAO

export 
Foil, Wing, HalfWing, projected_area, span, mean_aerodynamic_chord, 
Panel3D, Uniform, print_info, plot_setup, panel_coords, read_foil, 
mesh_horseshoes, mesh_wing, mesh_cambers, 
solve_case, streamlines,
plot_panels, plot_streamlines

#----------------------IMPORTS--------------------------------#

include("MathTools.jl")
include("FoilParametrization.jl")
include("LiftingLine.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum, fwddiff, fwddiv, tuparray, vectarray, tupvector, dot, linspace, cosine_dist
using .FoilParametrization: read_foil, cosine_foil, foil_camthick
using .LiftingLine
using StaticArrays
using LinearAlgebra
using Rotations


#-------------------------AIRCRAFT GUFF---------------------#

abstract type Aircraft end

#-------------------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

"""
Scales the coordinates of a Foil, usually to some chord length.
"""
scale_foil(foil :: Foil, chord) = chord * foil.coords

"""
Translates the coordinates of a Foil by (x, y, z).
"""
shift_foil(foil :: Foil, x, y, z) = [ x y z ] .+ foil.coords

"""
Returns a Foil with cosine spacing for a given number of points. 
"""
cut_foil(foil :: Foil, num) = Foil(cosine_foil(foil.coords, n = num))

"""
Computes the camber-thickness distribution of a Foil with cosine spacing..
"""
camber_thickness(foil :: Foil, num) = Foil(foil_camthick(cosine_foil(foil.coords), num + 1))

"""
Projects a Foil onto the x-z plane at y = 0.
"""
coordinates(foil :: Foil) = [ foil.coords[:,1] (zeros ∘ length)(foil.coords[:,1]) foil.coords[:,2] ]

#-----------------WING---------------------#

"""
Definition for a half-wing consisting of airfoils, span lengths, dihedrals, and sweep angles.
"""
struct HalfWing <: Aircraft
    foils :: Array{Foil} # Airfoil profiles
    chords :: Array{Float64} # Chord lengths (m)
    spans :: Array{Float64}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{Float64} # Dihedral angles (deg)
    sweeps :: Array{Float64} # Leading-edge sweep angles (deg)
    twists :: Array{Float64} # Twist angles (deg)
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
    taper_ratios = fwddiv(wing.chords)
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord.(wing.chords[1:end-1], taper_ratios)
    sum(macs .* areas) / sum(areas)
end

"""
Computes the weighted average (μ) of two vectors. 
"""
vector_point(x1, x2, μ) = x1 .+ μ .* (x2 .- x1)

"""
Divides two vectors into `n' sections with cosine spacing. TODO: Upgrade to generic spacing functions.
"""
chop_sections(set1, set2, n) = [ vector_point.(set1, set2, μ) for μ in cosine_dist(0.5, 1, n + 1) ][1:end-1]

"""
Divides a set of directional vectors into `n' sections with cosine spacing.
"""
coords_chopper(coords, n) = [ (chop_sections.(coords[1:end-1], coords[2:end], n)...)..., coords[end] ]

"""
Computes the leading edge coordinates of a HalfWing, with an option to flip the signs of the y-coordinates.
"""
function lead_wing(wing :: HalfWing, flip :: Bool = false)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

    sweeped_spans = [ 0; cumsum(spans .* tan.(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(spans .* tan.(dihedrals)) ]
    cum_spans = [ 0; cumsum(spans) ]
    
    SVector.(sweeped_spans, flip ? -cum_spans : cum_spans, dihedraled_spans)
end

"""
Computes the leading and trailing edge coordinates of a HalfWing, with an option to flip the signs of the y-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip :: Bool = false)
    chords = wing.chords
    twisted_chords = chords .* sin.(wing.twists)
    
    leading = lead_wing(wing, flip)
    trailing = SVector.(chords, (zeros ∘ length)(chords), twisted_chords) .+ leading

    leading, trailing
end

"""
Computes the coordinates of a HalfWing consisting of Foils and relevant geometric quantities. Requires specification of number of spanwise and chordwise panels. Optionally flips the signs of the y-coordinates.
"""
function wing_coords(wing :: HalfWing, chord_num :: Integer, span_num :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ cut_foil).(wing.foils, chord_num)
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num)
end

"""
Computes the coordinates of a HalfWing consisting of camber distributions of Foils and relevant geometric quantities. Requires specification of number of spanwise and chordwise panels. Optionally flips the signs of the y-coordinates.
"""
function camber_coords(wing :: HalfWing, chord_num :: Integer, span_num :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ camber_thickness).(wing.foils, chord_num)
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num)
end

"""
Zips leading and trailing edge coordinates.
"""
chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) in zip(lead, trail) ]

"""
Returns an array of cosine distributed coordinates 
"""
chord_chopper(coords, divs) = [ [ vector_point(chord[1,:], chord[2,:], μ) for μ in cosine_dist(0.5, 1, divs + 1) ] for chord in coords ]

"""
Applies spanwise spacing divisions on leading and trailing edge coordinates.
"""
span_chopper(lead, trail, div) = coords_chopper(lead, div), coords_chopper(trail, div)

"""
Chops the bounds of a wing (leading and trailing edge coordinates) into a wingbox for a given number of spanwise and chordwise panels.
"""
wing_chopper(lead, trail, span_num, chord_num) = chord_chopper(chord_sections(span_chopper(lead, trail, span_num)...), chord_num)

"""
Useless for now, but looks cool.
"""
panel(root_lead, root_trail, tip_trail, tip_lead) = [ root_lead  tip_lead;
                                                      root_trail tip_trail ]

"""
Converts an array of "wing-ordered" coordinates into panels.
"""
function make_panels(coords) 
    spanlist = vectarray.(coords)    
    secs1secs2 = zip(spanlist, spanlist[2:end])
    hcat([ Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) in secs1secs2 ]...)
end

"""
Meshes a HalfWing into panels meant for lifting-line analyses using horseshoe elements.
"""
mesh_horseshoes(obj :: HalfWing; span_num :: Integer = 1, chord_num :: Integer = 1, flip = false) = (make_panels ∘ wing_chopper)(wing_bounds(obj, flip)..., span_num, chord_num)

"""
Meshes a HalfWing into panels meant for full 3D analyses using doublet-source elements. Also works as a surface mesh. TODO: Tip meshing.
"""
mesh_wing(wing :: HalfWing; chord_num :: Integer = 1, span_num :: Integer = 1, flip = false) = (make_panels ∘ wing_coords)(wing, chord_num, span_num, flip)

"""
Meshes a HalfWing into panels meant for vortex lattice analyses using vortex rings.
"""
mesh_cambers(wing :: HalfWing; chord_num :: Integer = 1, span_num :: Integer = 1, flip = false) = (make_panels ∘ camber_coords)(wing, chord_num, span_num, flip)

"""
Performs an affine transformation on a list of coordinates.
"""
transform(coords, rotation, translation) = coords * rotation' .+ translation

"""
Performs an affine transformation on the leading and trailing edges of a HalfWing.
"""
transform(lead, trail; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = transform(lead, rotation, translation), transform(trail, rotation', translation)

"""
A composite type consisting of two HalfWings. `left' is meant to have its y-coordinates flipped in Cartesian representation.
"""
struct Wing <: Aircraft
    left :: HalfWing
    right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Union{Wing, HalfWing}) = aspect_ratio(span(wing), projected_area(wing))

"""
Computes the wing bounds for a Wing.
"""
function wing_bounds(wing :: Wing)
    left_lead, left_trail = wing_bounds(wing.left, true)
    right_lead, right_trail = wing_bounds(wing.right)

    leading = [ left_lead[end:-1:2,:]; right_lead ]
    trailing = [ left_trail[end:-1:2,:]; right_trail ]

    leading, trailing
end

"""
Meshes a Wing into panels meant for lifting-line analyses using horseshoe elements.
"""
function mesh_horseshoes(wing :: Wing; chord_num = 20, span_num = 3)
    left_panels = mesh_horseshoes(wing.left, chord_num = chord_num, span_num = span_num, flip = true)
    right_panels = mesh_horseshoes(wing.right, chord_num = chord_num, span_num = span_num)

    [ left_panels; right_panels ] 
end

"""
Meshes a HalfWing into panels meant for full 3D analyses using doublet-source elements. Also works as a surface mesh. TODO: Tip meshing.
"""
function mesh_wing(wing :: Wing; chord_num = 20, span_num = 3)
    left_panels = mesh_wing(wing.left, chord_num = chord_num, span_num = span_num, flip = true)
    right_panels = mesh_wing(wing.right, chord_num = chord_num, span_num = span_num)

    [ left_panels; right_panels ] 
end

"""
Meshes a HalfWing into panels meant for vortex lattice analyses using vortex rings.
"""
function mesh_cambers(wing :: Wing; chord_num = 20, span_num = 3)
    left_panels = mesh_cambers(wing.left, chord_num = chord_num, span_num = span_num, flip = true)
    right_panels = mesh_cambers(wing.right, chord_num = chord_num, span_num = span_num)

    [ left_panels; right_panels ] 
end

"""
Prints the relevant geometric characteristics of a HalfWing or Wing.
"""
function print_info(wing :: Union{Wing, HalfWing})
    println("Span: ", span(wing), " m")
    println("Area: ", projected_area(wing), " m²")
    println("MAC: ", mean_aerodynamic_chord(wing), " m")
    println("Aspect Ratio: ", aspect_ratio(wing))
end

function solve_case(wing :: Aircraft, uniform :: Uniform3D, r_ref = (0.25, 0, 0), ρ = 1.225; span_num = 15, chord_num = 5)

    # Make panels
    horseshoe_panels = mesh_horseshoes(wing, span_num = span_num, chord_num = chord_num)
    
    # Solve system
    Γs, horseshoes = solve_horseshoes(horseshoe_panels, uniform)

    # Compute forces
    geom_forces, geom_moments = dynamic_computations(Γs, horseshoes, uniform, r_ref, ρ)
    
    # Stability axes transformation
    force, moment = sum(geom_forces), sum(geom_moments)
    stable_forces, stable_moments = stability_axes(force, moment, uniform)
    drag = nearfield_drag(force, uniform)

    # Non-dimensionalisation parameters
    V = uniform.mag
    S = projected_area(wing)
    b = span(wing)
    c = mean_aerodynamic_chord(wing);

    # Print data
    print_dynamics(force, moment, drag, V, S, b, c, ρ)

    horseshoe_panels, horseshoes, Γs
end

plot_panels(panels :: Array{Panel3D}) = (tuparray ∘ panel_coords).(panels)
plot_streamlines(streams) = tupvector(streams)



end