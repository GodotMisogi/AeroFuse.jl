module AeroMDAO

export make_panels, sections, coordinates, print_info, plot_setup, Wing, HalfWing, Foil, Panel3D, horseshoe_collocation, horseshoe_point, horseshoe_vortex, Uniform, solve_case, projected_area, horseshoe_lines, transform

include("MathTools.jl")
# include("DoubletSource.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum, fwddiff, tuparray, dot, linspace
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

function wing_coords(wing :: HalfWing)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps
    root_chord, tip_chords = peel(wing.chords)
    root_twist, tip_twists = peel(wing.twists)

    sweeped_spans, dihedraled_spans, cum_spans = cumsum(spans .* tan.(sweeps)), cumsum(spans .* tan.(dihedrals)), cumsum(spans)
    twisted_chords = tip_chords .* sin.(tip_twists)

    leading_xyz = [ sweeped_spans cum_spans dihedraled_spans ]
    trailing_xyz = [ (tip_chords .+ sweeped_spans) (cum_spans) (dihedraled_spans .+ twisted_chords) ]

    leading = [ 0 0 0; leading_xyz ]
    trailing = [ root_chord 0 root_chord * sin(root_twist); trailing_xyz ]

    leading, trailing
end

chopper(xs, divs) = [ ([ [ x1 + μ*(x2 - x1) for μ in 0:1/divs:1 ][1:end-1] for (x2, x1) in zip(xs[2:end], xs) ]...); xs[end] ]

sections(lead, trail) = [ [ le'; te' ] for (le, te) in zip(eachrow(lead), eachrow(trail)) ]

coords_chopper(coords, divs) = [ chopper(coords[:,1], divs) chopper(coords[:,2], divs) chopper(coords[:,3], divs) ]

chord_chopper(lead, trail, chordwise_panels) = coords_chopper.(sections(lead, trail), chordwise_panels)

span_chopper(lead, trail, spanwise_panels) = coords_chopper(lead, spanwise_panels), coords_chopper(trail, spanwise_panels)

wing_chopper(lead, trail, spanwise_panels, chordwise_panels) = chord_chopper(coords_chopper(lead, spanwise_panels), coords_chopper(trail, spanwise_panels), chordwise_panels)

transform(lead, trail; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = lead * rotation' .+ translation, trail * rotation' .+ translation


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

function wing_coords(wing :: Wing)
    left_lead, left_trail = wing_coords(wing.left)
    right_lead, right_trail = wing_coords(wing.right)

    yflip!(left_lead)
    yflip!(left_trail)

    leading = [ left_lead[end:-1:2,:]; right_lead ]
    trailing = [ left_trail[end:-1:2,:]; right_trail ]

    leading, trailing
end

make_panels(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ], spanwise_panels = 1) = make_panels(wing_chopper(transform(wing_coords(obj)..., rotation = rotation, translation = translation)..., spanwise_panels)...)

sections(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = sections(transform(wing_coords(obj)..., rotation = rotation, translation = translation)...)

coordinates(obj :: Union{Wing, HalfWing}; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = coordinates(transform(wing_chopper(wing_coords(obj)...)..., rotation = rotation, translation = translation)...)


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
    α :: Float64 
    β :: Float64
end

Uniform(mag, α_deg, β_deg) = Uniform3D(mag, deg2rad(α_deg), deg2rad(β_deg))

freestream_to_cartesian(r, θ, φ) = r .* (cos(θ) * cos(φ), cos(θ) * sin(φ), sin(θ))

velocity(uni :: Uniform3D) = freestream_to_cartesian(uni.mag, uni.α, uni.β)

# Axes definitions
stability_axes(obj :: Union{Wing, HalfWing}, uni :: Uniform3D) = coordinates(obj, rotation = RotY{Float64}(-uni.α))
wind_axes(obj :: Union{Wing, HalfWing}, uni :: Uniform3D) = coordinates(obj, rotation = RotZY{Float64}(-uni.β, -uni.α))

abstract type Panel end

struct Panel3D <: Panel
    p1 :: SVector{3,Float64}
    p2 :: SVector{3,Float64}
    p3 :: SVector{3,Float64}
    p4 :: SVector{3,Float64}
end

collocation_point(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) / 4
panel_normal(panel :: Panel3D) = let p21 = panel.p2 .- panel.p1, p41 = panel.p4 .- panel.p1, p21_x_p41 = p21 × p41;
                                 p21_x_p41 / norm(p21_x_p41) end


weighted(x1, x2, wx) = (1 - wx) * x1 + wx * x2
weighted_point((x1, y1, z1), (x2, y2, z2), wx, wy, wz) = (weighted(x1, x2, wx), weighted(y1, y2, wy), weighted(z1, z2, wz))
quarter_point(p1, p2) = weighted_point(p1, p2, 1/4, 0, 1/2)
three_quarter_point(p1, p2) = weighted_point(p1, p2, 3/4, 0, 1/2)

horseshoe_line(p1, p2, p3, p4) = [ quarter_point(p1, p2); quarter_point(p4, p3) ]
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) .+ three_quarter_point(p4, p3) ) ./ 2

struct Line
    r1 :: SVector{3,Float64}
    r2 :: SVector{3,Float64}
end

function velocity(r, line :: Line, Γ, ε = 1e-6)
    r1 = r .- line.r1
    r2 = r .- line.r2
    r0 = line.r2 .- line.r1
    r1_x_r2 = r1 × r2

    any(<(ε), norm.([r1, r2, r1_x_r2])) ? [0,0,0] : Γ / (4π) * dot(r0, r1 / norm(r1) .- r2 / norm(r2)) * r1_x_r2 / norm(r1_x_r2) 
end


struct Horseshoe
    vortex_lines :: Array{Line}
end

function Horseshoe(vortex_line :: Line, uniform :: Uniform3D, bound :: Float64 = 30.)
    r1, r2 = vortex_line.r1, vortex_line.r2
    vel = bound .* velocity(uniform) ./ uniform.mag
    tline_1 = Line(r1 .+ vel, r1)
    tline_2 = Line(r2, r2 .+ vel)

    Horseshoe([ tline_1, vortex_line, tline_2 ])
end

horseshoe_vortex(panel :: Panel3D) = horseshoe_line(panel.p1, panel.p2, panel.p3, panel.p4)

horseshoe_collocation(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

horseshoe_lines(panel :: Panel3D, uniform :: Uniform3D) = Horseshoe(Line(horseshoe_vortex(panel)...), uniform)

velocity(r, horseshoe :: Horseshoe, Γ) = sum([ velocity(r, line, Γ) for line in horseshoe.vortex_lines ])

downwash_velocity(r, horseshoe :: Horseshoe, Γ) = let vels = [ velocity(r, line, Γ) for line in horseshoe.vortex_lines ]; first(vels) + last(vels) end

#---------------------------------Matrix setup--------------------------------------#

influence_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal) = dot(velocity(collocation_point, horseshoe, 1.), panel_normal)

induced_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal) = dot(downwash_velocity(collocation_point, horseshoe, 1.), panel_normal)

influence_matrix(colpoints, horseshoes :: Array{<: Horseshoe}, normals) = [ influence_coefficient(colpoint_i, horsey_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsey_j in horseshoes ]

induced_matrix(colpoints, horseshoes :: Array{<: Horseshoe}, normals)  = [ induced_coefficient(colpoint_i, horsey_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsey_j in horseshoes ]

boundary_condition(normals, velocity) = - [ dot(velocity, normal) for normal in normals ]

#-------------------------Force evaluations------------------------------------#

lift(Γ, Δy, speed, ρ = 1.225) = ρ * speed * Γ * Δy
induced_drag(Γ, Δy, w_ind, ρ = 1.225) = -ρ * w_ind * Γ * Δy

# Trefftz plane evaluations
# downwash(xi, xj, zi, zj) = -1/(2π) * (xj - xi) / ( (zj - zi)^2 + (xj - xi)^2 )
# trefftz_plane(horseshoes :: Array{Horseshoes}) = [ [ horseshoe.vortex_lines[1].r1 for horseshoe in horseshoes ]; 
#                                                    [ horseshoe.vortex_lines[3].r2 for horseshoe in horseshoes ] ]
# function downwash(horseshoes :: Array{Horseshoes}) 
#     coords = trefftz_plane(horseshoes)
#     carts = product(coords, coords)
#     [ [ i == j ? 0 : downwash(x1, x2, z1, z2) for (i, (x1,y1,z1)) in enumerate(coords) ] for (j, (x2,y2,z2)) in enumerate(coords) ] 
# end

function solve_case(panels :: Array{<: Panel}, uniform :: Uniform3D, ρ = 1.225) 
    
    horseshoes = [ horseshoe_lines(panel, uniform) for panel in panels ]
    colpoints = horseshoe_collocation.(panels)
    normals = panel_normal.(panels)
    horsies = horseshoe_vortex.(panels) 
    vel = velocity(uniform)

    Γs = influence_matrix(colpoints, horseshoes, normals) \ boundary_condition(normals, vel)
    w_inds = induced_matrix(colpoints, horseshoes, normals) * Γs
    Δys = (abs ∘ norm).([ line[2] .- line[1] for line in horsies ])

    Lift = sum(lift.(Γs, Δys, uniform.mag, ρ))
    Induced_drag = sum(induced_drag.(Γs, Δys, w_inds, ρ))
    Trefftz_drag = 

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