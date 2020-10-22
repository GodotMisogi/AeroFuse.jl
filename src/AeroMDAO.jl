module AeroMDAO

export make_panels, sections, coordinates, print_info, plot_setup, Wing, HalfWing, Foil, Panel3D, horseshoe_collocation, horseshoe_point, horseshoe_vortex, Uniform3D, solve_case

include("MathTools.jl")
# include("DoubletSource.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum, fwddiff, tuparray, dot
# using .DoubletSource: Point2D, Point3D, Panel2D, collocation_point
using StaticArrays
using LinearAlgebra
using Rotations


yflip!(xs) = xs[:,2] .= -xs[:,2]
plot_setup(coords) = tuparray(coords)

abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

coordinates(foil :: Foil; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = foil.coords * rotation' .+ translation

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
    taper_ratios = wing.chords[2:end] ./ wing.chords[1:end-1]
    areas = mean_chords .* wing.spans
    macs = mean_aerodynamic_chord.(wing.chords[1:end-1], taper_ratios)
    sum(macs .* areas) / sum(areas)
end

function wing_coords(wing :: HalfWing, rotation, translation)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps
    root_chord, tip_chords = peel(wing.chords)
    root_twist, tip_twists = peel(wing.twists)

    sweeped_spans, dihedraled_spans, cum_spans = cumsum(spans .* sin.(sweeps)), cumsum(spans .* sin.(dihedrals)), cumsum(spans)
    twisted_chords = tip_chords .* sin.(tip_twists)

    leading_xyz = [ sweeped_spans cum_spans dihedraled_spans ]
    trailing_xyz = [ (tip_chords .+ sweeped_spans) (cum_spans) (dihedraled_spans .+ twisted_chords) ]

    leading = [ 0 0 0; leading_xyz ] * rotation' .+ translation
    trailing = [ root_chord 0 root_chord * sin(root_twist); trailing_xyz ] * rotation' .+ translation

    leading, trailing
end

sections(lead, trail) = [ [ le'; te' ] for (le, te) in zip(eachrow(lead), eachrow(trail)) ]
coordinates(lead, trail) = [ lead[end:-1:1,:]; trail ]

function make_panels(lead, trail)
    lead, trail = tuparray(lead), tuparray(trail)
    
    panel_coords = zip(lead, trail[1:end-1], trail[2:end], lead[2:end])

    [ Panel3D(x...) for x in panel_coords ]
end

struct Wing <: Aircraft
    left :: HalfWing
    right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Union{Wing, HalfWing}) = aspect_ratio(span(wing), projected_area(wing))

function wing_coords(wing :: Wing, rotation, translation)
    left_lead, left_trail = wing_coords(wing.left, rotation, translation)
    right_lead, right_trail = wing_coords(wing.right, rotation, translation)

    yflip!(left_lead)
    yflip!(left_trail)

    leading = [ left_lead[end:-1:2,:]; right_lead ]
    trailing = [ left_trail[end:-1:2,:]; right_trail ]

    leading, trailing
end

make_panels(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = make_panels(wing_coords(obj, rotation, translation)...)
sections(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = sections(wing_coords(obj, rotation, translation)...)
coordinates(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = coordinates(wing_coords(obj, rotation, translation)...)

function print_info(wing :: Union{Wing, HalfWing})
    println("Span: ", span(wing), " m")
    println("Area: ", projected_area(wing), " m²")
    println("MAC: ", mean_aerodynamic_chord(wing), " m")
    println("Aspect Ratio: ", aspect_ratio(wing))
end

#--------------Lifting line code-------------------#

abstract type Laplace end

struct Uniform3D <: Laplace
    mag :: Float64
    alpha :: Float64 
    beta :: Float64
end

velocity(uni :: Uniform3D) = let alpha = deg2rad(uni.alpha), beta = deg2rad(uni.beta); 
    uni.mag .* (cos(alpha) * cos(beta), cos(alpha) * sin(beta), sin(alpha)) end

struct Doublet2D <: Laplace
    str :: Float64
    x0 :: Float64
    y0 :: Float64 
end

potential(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)

abstract type Panel end

struct Panel3D <: Panel
    p1 :: SVector{3,Float64}
    p2 :: SVector{3,Float64}
    p3 :: SVector{3,Float64}
    p4 :: SVector{3,Float64}
end

collocation_point(panel :: Panel3D) = [ panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4 ] / 4
panel_normal(panel :: Panel3D) = let p21 = panel.p2 .- panel.p1, p41 = panel.p4 .- panel.p1, p21_x_p41 = p21 × p41;
                                 p21_x_p41 / norm(p21_x_p41) end
quarter_point((x1, y1, z1), (x2, y2, z2)) = ( (3 * x1 + x2) / 4, y1, (z1 + z2) / 2 )


horseshoe_line(p1, p2, p3, p4) = [ quarter_point(p1, p2); quarter_point(p4, p3) ]
collocation_point(p1, p2, p3, p4) = ( ( 3 * p1[1] + p2[1] + p3[1] + 3 * p4[1]) / 8, (p1[2] + p3[2]) / 2, (p1[3] + p2[3] + p3[3] + p4[3]) / 4 )

horseshoe_collocation(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

struct Line
    r1 :: SVector{3,Float64}
    r2 :: SVector{3,Float64}
end

horseshoe_vortex(panel :: Panel3D) = horseshoe_line(panel.p1, panel.p2, panel.p3, panel.p4)

velocity(line :: Line, r, Γ, ε = 1e-6) = let r1 = r .- line.r1, r2 = r .- line.r2, r1_x_r2 = r1 × r2;
    any(<(ε), norm.([r1, r2, r1_x_r2])) ? [0,0,0] : Γ/(4π) * r1_x_r2 / norm(r1_x_r2) .* r1 .- r2 .* (r1 / norm(r1) .- r2 / norm(r2)) end

mutable struct Horseshoe
    vortex_lines :: Array{Line}
    function Horseshoe(vortex_line :: Line, angle :: Float64, bound :: Float64 = 1e5)
        r1, r2 = vortex_line.r1, vortex_line.r2
        tline_1 = Line([ bound, r1[2], r1[3] * (sin ∘ deg2rad)(angle) ], r1)
        tline_2 = Line(r2, [ bound, r2[2], r2[3] * (sin ∘ deg2rad)(angle) ])

        new([vortex_line, tline_1, tline_2])
    end
end

horseshoe_lines(panel :: Panel3D, angle = 0.) = Horseshoe(Line(horseshoe_vortex(panel)...), angle)

velocity(horseshoe :: Horseshoe, r, Γ) = sum(velocity.(horseshoe.vortex_lines, r, Γ))
downwash_velocity(horseshoe :: Horseshoe, r, Γ) = let vels = velocity.(horseshoe.vortex_lines, r, Γ, 1); first(vels) + last(vels) end

#---------------------------------Matrix setup--------------------------------------#

influence_coefficient(horseshoe :: Horseshoe, panel_2 :: Panel3D) = dot(velocity(horseshoe, horseshoe_collocation(panel_2), 1.), panel_normal(panel_2))

induced_coefficient(horseshoe :: Horseshoe, panel_2 :: Panel3D) = dot(downwash_velocity(horseshoe, horseshoe_collocation(panel_2), 1.), panel_normal(panel_2))

influence_matrix(panels :: Array{<: Panel}) = [ influence_coefficient(horsey_j, panel_i) for horsey_j in horseshoe_lines.(panels), panel_i in panels ]

induced_matrix(panels :: Array{<: Panel}) = [ induced_coefficient(horsey_j, panel_i) for horsey_j in horseshoe_lines.(panels), panel_i in panels ]

boundary_condition(panels :: Array{<: Panel}, uniform :: Uniform3D) = - [ dot(velocity(uniform), panel_normal(panel)) for panel in panels ]

#-------------------------Force evaluations------------------------------------#

lift(Γ, Δy, speed, ρ = 1.225) = ρ * speed * Γ * Δy
induced_drag(Γ, Δy, w_ind, speed, ρ = 1.225) = -ρ * w_ind * Γ * Δy

function solve_case(panels :: Array{<: Panel}, uniform :: Uniform3D) 
    
    Γs = influence_matrix(panels) \ boundary_condition(panels, uniform)
    
    w_ind = induced_matrix(panels) * Γs

    horsies = horseshoe_vortex.(panels)
    Δys = (abs ∘ norm).([ line[2] .- line[1] for line in horsies ])
    gammyind = zip(Γs, Δys, w_ind)
    println(collect(gammyind))
    Lift = sum([ lift(Γ, Δy, uniform.mag) for (Γ, Δy, w_ind) in gammyind ] )
    Induced_drag = sum([ induced_drag(pt..., uniform.mag) for pt in gammyind ])

    Lift, Induced_drag
end

end


# Kleisli-ish composition to transport global coordinates?
# PlanePos = Pair(Aircraft :: Aircraft, Point3D)

# function mesh_wing(wing :: HalfWing, span_panels = 8)
#     lead, trail = wing_coords(wing)

#     foils = wing.foils
#     chords = wing.chords
#     spans = wing.spans

#     span_increments = (fwddiff ∘ cumsum)(spans) / span_panels


#     coords = chords .* (:coords .<< foils)


# end