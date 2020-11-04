module AeroMDAO

export Foil, Wing, HalfWing, Panel3D, Uniform, coordinates, print_info, plot_setup, transform, cut_foil, panel_coords, read_foil, mesh_horseshoes, mesh_wing, mesh_cambers, horseshoe_collocation, horseshoe_point, horseshoe_vortex, Uniform, solve_case, projected_area, horseshoe_lines

include("MathTools.jl")
include("FoilParametrization.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum, fwddiff, fwddiv, tuparray, vectarray, dot, linspace
using .FoilParametrization: read_foil, cosine_foil, foil_camthick
using StaticArrays
using LinearAlgebra
using Rotations

plot_setup(coords) = tuparray(coords)

abstract type Aircraft end

#----------------AIRFOIL----------------------#

"""
Airfoil structure consisting of foil coordinates as an array of points.
"""
struct Foil <: Aircraft
    coords :: Array{<: Real, 2} # The foil profile as an array of coordinates, must be in Selig format.
end

scale_foil(foil :: Foil, chord) = chord * foil.coords
shift_foil(foil :: Foil, x, y, z) = [ x y z ] .+ foil.coords 
cut_foil(foil :: Foil, num) = Foil(cosine_foil(foil.coords, n = num))
cambers(foil :: Foil, num) = Foil(foil_camthick(cosine_foil(foil.coords), num + 1))

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

vector_point(x1, x2, μ) = x1 .+ μ .* (x2 .- x1)
chop_sections(set1, set2, divs) = [ vector_point.(set1, set2, μ) for μ in 0:1/divs:1 ][1:end-1]

function lead_wing(wing :: HalfWing, flip :: Bool = false)
    spans, dihedrals, sweeps = wing.spans, wing.dihedrals, wing.sweeps

    sweeped_spans = [ 0; cumsum(spans .* tan.(sweeps)) ]
    dihedraled_spans = [ 0; cumsum(spans .* tan.(dihedrals)) ]
    cum_spans = [ 0; cumsum(spans) ]
    
    SVector.(sweeped_spans, flip ? -cum_spans : cum_spans, dihedraled_spans)
end

function wing_bounds(wing :: HalfWing, flip :: Bool = false)
    chords = wing.chords
    twisted_chords = chords .* sin.(wing.twists)
    
    leading = lead_wing(wing, flip)
    trailing = SVector.(chords, (zeros ∘ length)(chords), twisted_chords) .+ leading

    leading, trailing
end

coords_chopper(coords, divs) = [ (chop_sections.(coords[1:end-1], coords[2:end], divs)...)..., coords[end] ]

function wing_coords(wing :: HalfWing, chordwise_panels :: Integer, spanwise_panels :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* coordinates.(cut_foil.(wing.foils, chordwise_panels))
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, spanwise_panels)
end

function camber_coords(wing :: HalfWing, chordwise_panels :: Integer, spanwise_panels :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ cambers).(wing.foils, chordwise_panels)
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, spanwise_panels)
end

chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) in zip(lead, trail) ]

chord_chopper(coords, divs) = [ [ vector_point(chord[1,:], chord[2,:], μ) for μ in 0:1/divs:1 ] for chord in coords ]

span_chopper(lead, trail, div) = coords_chopper(lead, div), coords_chopper(trail, div)

wing_chopper(lead, trail, spanwise_panels, chordwise_panels) = chord_chopper(chord_sections(span_chopper(lead, trail, spanwise_panels)...), chordwise_panels)

coordinates(lead, trail) = [ lead[end:-1:1,:]; trail ]

panel(root_lead, root_trail, tip_trail, tip_lead) = [ root_lead  tip_lead;
                                                      root_trail tip_trail ]

function make_panels(coords) 
    spanlist = vectarray.(coords)    
    secs1secs2 = zip(spanlist, spanlist[2:end])
    hcat([ Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) in secs1secs2 ]...)
end

mesh_horseshoes(obj :: HalfWing; spanwise_panels = 1, chordwise_panels = 1, flip = false) = (make_panels ∘ wing_chopper)(wing_bounds(obj, flip)..., spanwise_panels, chordwise_panels)

mesh_wing(wing :: HalfWing; chordwise_panels :: Integer = 20, spanwise_panels:: Integer = 3, flip = false) = (make_panels ∘ wing_coords)(wing, chordwise_panels, spanwise_panels, flip)

mesh_cambers(wing :: HalfWing; chordwise_panels :: Integer = 20, spanwise_panels:: Integer = 3, flip = false) = (make_panels ∘ camber_coords)(wing, chordwise_panels, spanwise_panels, flip)

transform(lead, trail; rotation = one(RotMatrix{3, Float64}), translation = [ 0 0 0 ]) = lead * rotation' .+ translation, trail * rotation' .+ translation

struct Wing <: Aircraft
    left :: HalfWing
    right :: HalfWing
end

span(wing :: Wing) = span(wing.left) + span(wing.right)
projected_area(wing :: Wing) = projected_area(wing.left) + projected_area(wing.right)
mean_aerodynamic_chord(wing :: Wing) = (mean_aerodynamic_chord(wing.left) + mean_aerodynamic_chord(wing.right)) / 2
aspect_ratio(wing :: Union{Wing, HalfWing}) = aspect_ratio(span(wing), projected_area(wing))

function wing_bounds(wing :: Wing)
    left_lead, left_trail = wing_bounds(wing.left, true)
    right_lead, right_trail = wing_bounds(wing.right)

    leading = [ left_lead[end:-1:2,:]; right_lead ]
    trailing = [ left_trail[end:-1:2,:]; right_trail ]

    leading, trailing
end

function mesh_horseshoes(wing :: Wing; chordwise_panels = 20, spanwise_panels = 3)
    left_panels = mesh_horseshoes(wing.left, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels, flip = true)
    right_panels = mesh_horseshoes(wing.right, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels)

    [ left_panels; right_panels ] 
end

function mesh_wing(wing :: Wing; chordwise_panels = 20, spanwise_panels = 3)
    left_panels = mesh_wing(wing.left, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels, flip = true)
    right_panels = mesh_wing(wing.right, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels)

    [ left_panels; right_panels ] 
end

function mesh_cambers(wing :: Wing; chordwise_panels = 20, spanwise_panels = 3)
    left_panels = mesh_cambers(wing.left, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels, flip = true)
    right_panels = mesh_cambers(wing.right, chordwise_panels = chordwise_panels, spanwise_panels = spanwise_panels)

    [ left_panels; right_panels ] 
end

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
    #    p1 —→—— p4
    #    |      |
    #   ↓      ↓
    #  |      |
    # p2 —→— p3
end

panel_coords(panel :: Panel3D) = [ panel.p1'; panel.p2'; panel.p3'; panel.p4'; panel.p1' ] 

collocation_point(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) / 4
panel_normal(panel :: Panel3D) = let p21 = panel.p2 .- panel.p1, 
                                     p41 = panel.p4 .- panel.p1, 
                                     p21_x_p41 = p21 × p41;
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

struct VortexArray
    vortex_lines :: Array{Line}
end

function Horseshoe(vortex_line :: Line, uniform :: Uniform3D, bound :: Float64 = 30.)
    r1, r2 = vortex_line.r1, vortex_line.r2
    vel = bound .* velocity(uniform) ./ uniform.mag
    trailing_vortex_line_1 = Line(r1 .+ vel, r1)
    trailing_vortex_line_2 = Line(r2, r2 .+ vel)

    VortexArray([ trailing_vortex_line_1, vortex_line, trailing_vortex_line_2 ])
end

horseshoe_vortex(panel :: Panel3D) = horseshoe_line(panel.p1, panel.p2, panel.p3, panel.p4)

horseshoe_collocation(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

horseshoe_lines(panel :: Panel3D, uniform :: Uniform3D) = Horseshoe(Line(horseshoe_vortex(panel)...), uniform)

function vortex_ring(p1, p2, p3, p4)
    v1, v4 = SVector(quarter_point(p1, p2)...), SVector(quarter_point(p4, p3)...)
    v2, v3 = v1 .+ p2 .- p1, v4 .+ p3 .- p4
    [ v1, v2, v3, v4 ]
end 

vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

function VortexArray(panel :: Panel3D)
    v1, v2, v3, v4 = vortex_ring(panel)

    vortex_line = Line(v4, v1)
    left_line = Line(v1, v2)
    right_line = Line(v2, v3)
    back_line = Line(v3, v4)

    VortexArray([ left_line, vortex_line, back_line, right_line ])
end

velocity(r, vortex_line :: VortexArray, Γ) = sum([ velocity(r, line, Γ) for line in vortex_line.vortex_lines ])

downwash_velocity(r, vortex_line :: VortexArray, Γ) = let vels = [ velocity(r, line, Γ) for line in vortex_line.vortex_lines ]; first(vels) + last(vels) end

#---------------------------------Matrix setup--------------------------------------#

influence_coefficient(collocation_point, vortex_line :: VortexArray, panel_normal) = dot(velocity(collocation_point, vortex_line, 1.), panel_normal)

induced_coefficient(collocation_point, vortex_line :: VortexArray, panel_normal) = dot(downwash_velocity(collocation_point, vortex_line, 1.), panel_normal)

influence_matrix(colpoints, vortex_lines :: Array{<: VortexArray}, normals) = [ influence_coefficient(colpoint_i, horsey_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsey_j in vortex_lines ]

induced_matrix(colpoints, vortex_lines :: Array{<: VortexArray}, normals)  = [ induced_coefficient(colpoint_i, horsey_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsey_j in vortex_lines ]

boundary_condition(normals, velocity) = - [ dot(velocity, normal) for normal in normals ]

function wake_panel(panel :: Panel3D, uniform :: Uniform3D, bound = 10.)
    wake = bound * norm(panel.p2 .- panel.p1) .* velocity(uniform) ./ uniform.mag
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

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

#------------------------Case setup and solution--------------------------#

function accumap(f, n, xs)
    data = [ xs ]
    while n > 0
        println(n)
        ys = map(f, xs)
        data = [ data..., ys ]
        xs = ys
        n = n - 1
    end
    return hcat(data...)
end

function solve_case(panels :: Array{Panel3D}, uniform :: Uniform3D, wake_num = 10, ρ = 1.225) 
    
    vel = velocity(uniform)

    wake_panels = permutedims(accumap(panel -> wake_panel(panel, uniform), wake_num, panels[end,:]))
    panels = [ panels; 
               wake_panels ]

    horseshoes = [ horseshoe_lines(panel, uniform) for panel in panels ][:]
    horseshoes = VortexArray.(panels)[:]
    colpoints = horseshoe_collocation.(panels)[:]
    normals = panel_normal.(panels)[:]
    horsies = horseshoe_vortex.(panels)[:]

    inf = influence_matrix(colpoints, horseshoes, normals)
    boco = boundary_condition(normals, vel)

    Γs =  inf \ boco 

    w_inds = induced_matrix(colpoints, horseshoes, normals) * Γs
    Δys = (abs ∘ norm).([ line[2] .- line[1] for line in horsies ])

    Lifts = lift.(Γs, Δys, uniform.mag, ρ)
    Induced_drags = induced_drag.(Γs, Δys, w_inds, ρ)
    # Trefftz_drag = 

    Lifts, Induced_drags, panels
end


end