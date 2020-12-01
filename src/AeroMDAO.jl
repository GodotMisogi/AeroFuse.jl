module AeroMDAO

export 
Foil, Wing, HalfWing, make_wing, projected_area, span, mean_aerodynamic_chord, 
Panel3D, Freestream, info, print_info, plot_setup, panel_coords, panel_normal, read_foil, 
mesh_horseshoes, mesh_wing, mesh_cambers, lead_wing, wing_bounds,
solve_case, streamlines, tupvector, make_panels, vlmesh_wing, transform, paneller,
plot_panels, plot_streamlines, horseshoe_lines, aerodynamic_coefficients, print_dynamics,
pressure_coefficient, Uniform2D, make_2Dpanels, grid_data, split_panels, collocation_point

#----------------------IMPORTS--------------------------------#

include("MathTools.jl")
include("FoilParametrization.jl")
include("VortexLattice.jl")
include("DoubletSource.jl")
include("NonDimensionalization.jl")

using .MathTools: fwdsum, fwddiv, tuparray, vectarray, tupvector, linspace, cosine_dist, midgrad
using .FoilParametrization: read_foil, cosine_foil, foil_camthick
using .VortexLattice
using .DoubletSource
using .NonDimensionalization

import Base: *, +
using Base.Iterators: peel
using StaticArrays
using TimerOutputs
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
# chords :: Array{Float64} # Chord lengths (m)
struct HalfWing <: Aircraft
    foils :: Array{Foil} # Airfoil profiles
    chords :: Array{<: Real}
    spans :: Array{<: Real}  # Leading-edge to leading-edge distance between foils (m)
    dihedrals :: Array{<: Real} # Dihedral angles (deg)
    sweeps :: Array{<: Real} # Leading-edge sweep angles (deg)
    twists :: Array{<: Real} # Twist angles (deg)
    HalfWing(foils, chords, spans, dihedrals, sweeps, twists) = new(foils, chords, spans, deg2rad.(dihedrals), deg2rad.(sweeps), -deg2rad.(twists)) # Convert to radians
end

"""
Computes the aspect ratio given a span and area.
"""
aspect_ratio(span, area) = span^2 / area

"""
Computes the mean geometric chord given a span and area.
"""
mean_geometric_chord(span, area) = area / span

"""
Computes the taper ratio given root and tip chords.
"""
taper_ratio(root_chord, tip_chord) = tip_chord / root_chord

"""
Computes the area given a span and chord.
"""
area(span, chord) = span * chord

"""
Computes the mean aerodynamic chord, essentially a root-mean-square chord, given a root chord and taper ratio.
"""
mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)

"""
Computes the quarter-chord length given a chord.
"""
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
vector_point(x1, x2, μ) = x1 .+ μ * (x2 .- x1)

"""
Divides two vectors into `n` sections with cosine spacing. TODO: Upgrade to generic spacing functions.
"""
chop_sections(set1, set2, n) = [ vector_point.(set1, set2, μ) for μ ∈ cosine_dist(0.5, 1., n + 1) ][1:end-1]

"""
Divides a set of directional vectors into `n` sections with cosine spacing.
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
    
    leading = SVector.(sweeped_spans, flip ? -cum_spans : cum_spans, dihedraled_spans)

    flip ? leading[end:-1:1] : leading
end

"""
Computes the leading and trailing edge coordinates of a HalfWing, with an option to flip the signs of the y-coordinates.
"""
function wing_bounds(wing :: HalfWing, flip :: Bool = false)
    chords = wing.chords
    twisted_chords = chords .* sin.(wing.twists)
    
    leading = lead_wing(wing, flip)
    trailing = SVector.(chords, (zeros ∘ length)(chords), twisted_chords) 

    leading, (flip ? trailing[end:-1:1] : trailing) .+ leading
end

"""
Computes the coordinates of a HalfWing consisting of Foils and relevant geometric quantities. Requires specification of number of spanwise and chordwise panels. Optionally flips the signs of the y-coordinates.
"""
function wing_coords(wing :: HalfWing, span_num :: Integer, chord_num :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ cut_foil).(wing.foils, chord_num)
    
    scaled_foils = flip ? scaled_foils[end:-1:1] : scaled_foils
    twists = flip ? wing.twists[end:-1:1] : wing.twists

    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num)
end

"""
Computes the coordinates of a HalfWing consisting of camber distributions of Foils and relevant geometric quantities. Requires specification of number of spanwise and chordwise panels. Optionally flips the signs of the y-coordinates.
"""
function camber_coords(wing :: HalfWing, span_num :: Integer, chord_num :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ camber_thickness).(wing.foils, chord_num)

    scaled_foils = flip ? scaled_foils[end:-1:1] : scaled_foils
    twists = flip ? wing.twists[end:-1:1] : wing.twists

    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) ∈ zip(scaled_foils, twists, leading_xyz) ]

    coords_chopper(foil_coords, span_num)
end

"""
Zips leading and trailing edge coordinates into an array of arrays.
"""
chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) ∈ zip(lead, trail) ]

"""
Returns an array of cosine distributed coordinates 
"""
chord_chopper(coords, n) = [ [ vector_point(chord[1,:], chord[2,:], μ) for μ ∈cosine_dist(0.5, 1., n + 1)  ] for chord ∈ coords ]

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
    hcat((Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) ∈ secs1secs2 )...)
end

"""
Meshes a HalfWing into panels based on the chord lines of the airfoil sections, meant for lifting-line analyses using horseshoe elements.
"""
mesh_horseshoes(obj :: HalfWing, span_num :: Integer, chord_num :: Integer; flip = false) = (make_panels ∘ wing_chopper)(wing_bounds(obj, flip)..., span_num, chord_num)

"""
Meshes a HalfWing into panels meant for full 3D analyses using doublet-source elements. Also works as a surface mesh. TODO: Tip meshing.
"""
mesh_wing(wing :: HalfWing, span_num :: Integer, chord_num :: Integer; flip = false) = (make_panels ∘ wing_coords)(wing, span_num, chord_num, flip)

"""
Meshes a HalfWing into panels based on the camber distributions of the airfoil sections, meant for vortex lattice analyses.
"""
mesh_cambers(wing :: HalfWing, span_num :: Integer, chord_num :: Integer; flip = false) = (make_panels ∘ camber_coords)(wing, span_num, chord_num, flip)

# """
# Performs an affine transformation on a list of coordinates.
# """
# transform(coords, rotation, translation) = coords * rotation' .+ translation

# """
# Performs an affine transformation on the leading and trailing edges of a HalfWing.
# """
# transform(lead, trail; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = transform(lead, rotation', translation), transform(trail, rotation', translation)

"""
A composite type consisting of two HalfWings. `left' is meant to have its y-coordinates flipped ∈ Cartesian representation.
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
function mesh_horseshoes(wing :: Wing, span_num, chord_num)
    left_panels = mesh_horseshoes(wing.left, span_num, chord_num, flip = true)
    right_panels = mesh_horseshoes(wing.right, span_num, chord_num)

    [ left_panels right_panels ] 
end

"""
Meshes a Wing into panels meant for full 3D analyses using doublet-source elements. Also works as a surface mesh. TODO: Tip meshing.
"""
function mesh_wing(wing :: Wing, span_num, chord_num)
    left_panels = mesh_wing(wing.left, span_num, chord_num, flip = true)
    right_panels = mesh_wing(wing.right, span_num, chord_num)

    [ left_panels right_panels ] 
end

"""
Meshes a Wing into panels meant for vortex lattice analyses.
"""
function mesh_cambers(wing :: Wing, span_num, chord_num)
    left_panels = mesh_cambers(wing.left, span_num, chord_num, flip = true)
    right_panels = mesh_cambers(wing.right, span_num, chord_num)

    [ left_panels right_panels ] 
end

"""
Meshes a Half/Wing for vortex lattice analyses by generating horseshoe and camber panels.
"""
vlmesh_wing(wing :: Union{Wing, HalfWing}, span_num, chord_num) = mesh_horseshoes(wing, span_num, chord_num), mesh_cambers(wing, span_num, chord_num)

function paneller(wing :: Union{Wing, HalfWing}, span_num, chord_num; rotation = one(RotMatrix{3, Float64}), translation = SVector(0,0,0))
    horseshoes, cambers = vlmesh_wing(wing, span_num, chord_num) 
    [ transform(panel, rotation = rotation, translation = translation) for panel in horseshoes ],
    [ transform(panel, rotation = rotation, translation = translation) for panel in cambers ]
end

"""
Prints the relevant geometric characteristics of a HalfWing or Wing.
"""
function print_info(wing :: Union{Wing, HalfWing})
    data = info(wing)
    println("Span: ", data[1], " m")
    println("Area: ", data[2], " m²")
    println("MAC: ", data[3], " m")
    println("Aspect Ratio: ", data[4])
end

info(wing :: Union{Wing, HalfWing}) =
    span(wing), projected_area(wing), mean_aerodynamic_chord(wing), aspect_ratio(wing)

function solve_case(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, uniform :: Uniform3D, Ω = SVector(0., 0., 0.), r_ref = SVector(0.25, 0., 0.), ρ = 1.225; print = true, symmetry = false)
    vel = VortexLattice.velocity(uniform)

    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels, camber_panels, uniform, Ω, symmetry)

    # Compute near-field forces
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, uniform, Ω, r_ref, ρ)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, uniform.mag, VortexLattice.area.(camber_panels))
    # println(cps)

    force, moment = sum(geom_forces), sum(geom_moments)
    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, Ω, uniform)
    drag = nearfield_drag(force, uniform)

    # @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, uniform, r_ref, ρ)

    # # Compute non-dimensional coefficients
    # @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, uniform.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    force, drag, moment, horseshoes, Γs
    # horseshoe_panels, camber_panels, horseshoes, Γs
    # coeffs
end

function solve_case(wing :: Union{Wing, HalfWing}, uniform :: Uniform3D, Ω = SVector(0., 0., 0.), r_ref = SVector(0.25, 0., 0.), ρ = 1.225; span_num :: Integer = 15, chord_num :: Integer = 5, print = true)
    vel = VortexLattice.velocity(uniform)

    # Experimental: Symmetry condition
    symmetry = false
    # if typeof(wing) == Wing
    #     println("Symmetric Case")
    #     symmetry, wing = wing.left === wing.right ? (true, wing.right) : (false, wing)
    # end

    # reset_timer!()

    
    # Compute panels
    @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
    
    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels, camber_panels, uniform, Ω, symmetry)

    # Compute near-field forces
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, uniform, Ω, r_ref, ρ)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, uniform.mag, VortexLattice.area.(camber_panels))
    # println(cps)

    force, moment = sum(geom_forces), sum(geom_moments)
    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, Ω, uniform)
    drag = nearfield_drag(force, uniform)

    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, uniform, r_ref, ρ)

    # Compute non-dimensional coefficients
    @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, uniform.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, uniform.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    # # Optionally print data
    @timeit "Printing" if print 
        println("Nearfield:") 
        print_dynamics(nearfield_coeffs...)
        println("\nFarfield:")
        print_dynamics(farfield_coeffs...)
    end

    horseshoe_panels, camber_panels, horseshoes, Γs
    # nearfield_coeffs
end


# solve_case(wings :: Array{Aircraft}, uniform :: Uniform3D, r_ref = (0.25, 0, 0), ρ = 1.225; span_num = 15, chord_num = 5, print = true) =  

plot_panels(panels :: Array{Panel3D}) = (tupvector ∘ panel_coords).(panels)
plot_streamlines(streams) = tupvector(streams)

function solve_case(panels :: Array{Panel2D}, uniform :: Uniform2D)
    # Compute doublet and source strengths
    φs, σs = solve_strengths(panels, uniform)
    freestream_speed = norm(DoubletSource.velocity(uniform))
    pressure_coeffs = pressure_coefficient.(freestream_speed, panel_velocities(panels, uniform, φs[1:end-1]))

    # Make panels for plotting
    # dub_src_panels = DoubletSourcePanel2D.(panels, φs[1:end-1], σs, pressure_coeffs)

    # Compute lift coefficient
    diff_pans = [ panel_dist(panel_1, panel_2) for (panel_1, panel_2) ∈ (collect ∘ eachrow ∘ midgrad)(panels) ]
    cl = sum(lift_coefficient.(pressure_coeffs, diff_pans ./ 2, panel_angle.(panels)))

    cl
end

end