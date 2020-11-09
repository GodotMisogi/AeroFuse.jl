module AeroMDAO

export Foil, Wing, HalfWing, Panel3D, Uniform, coordinates, print_info, plot_setup, transform, cut_foil, panel_coords, read_foil, mesh_horseshoes, mesh_wing, mesh_cambers, vortex_collocation, horseshoe_point, bound_leg, Uniform, solve_horseshoes, solve_vortex_legs, projected_area, span, mean_aerodynamic_chord, horseshoe_lines, print_dynamics, streamlines

include("MathTools.jl")
include("FoilParametrization.jl")

import Base: *, +
using Base.Iterators: peel
using .MathTools: fwdsum, fwddiff, fwddiv, tuparray, vectarray, dot, linspace, cosine_dist
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
camber(foil :: Foil, num) = Foil(foil_camthick(cosine_foil(foil.coords), num + 1))

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
chop_sections(set1, set2, divs) = [ vector_point.(set1, set2, μ) for μ in cosine_dist(0.5, 1, divs + 1) ][1:end-1]

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
    
    scaled_foils = wing.chords .* (coordinates ∘ cut_foil).(wing.foils, chordwise_panels)
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, spanwise_panels)
end

function camber_coords(wing :: HalfWing, chordwise_panels :: Integer, spanwise_panels :: Integer, flip :: Bool = false)
    leading_xyz = lead_wing(wing, flip)
    
    scaled_foils = wing.chords .* (coordinates ∘ camber).(wing.foils, chordwise_panels)
    foil_coords = [ coords * RotY(-twist)' .+ section' for (coords, twist, section) in zip(scaled_foils, wing.twists, leading_xyz) ]

    coords_chopper(foil_coords, spanwise_panels)
end

chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) in zip(lead, trail) ]

chord_chopper(coords, divs) = [ [ vector_point(chord[1,:], chord[2,:], μ) for μ in cosine_dist(0.5, 1, divs + 1) ] for chord in coords ]

span_chopper(lead, trail, div) = coords_chopper(lead, div), coords_chopper(trail, div)

wing_chopper(lead, trail, spanwise_panels, chordwise_panels) = chord_chopper(chord_sections(span_chopper(lead, trail, spanwise_panels)...), chordwise_panels)

coordinates(lead, trail) = [ lead[end:-1:1,:]; trail ]

# Useless for now, but looks cool
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

freestream_to_cartesian(r, θ, φ) = r .* (cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

velocity(uni :: Uniform3D) = freestream_to_cartesian(uni.mag, uni.α, uni.β)


# Axes definitions
stability_axes(coords, uni :: Uniform3D) = RotY{Float64}(uni.α) * coords
wind_axes(coords, uni :: Uniform3D) = RotZY{Float64}(uni.β, uni.α) * coords

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

area(panel :: Panel3D) = (abs ∘ norm)((panel.p2 - panel.p1) .* (panel.p3 - panel.p2))
panel_coords(panel :: Panel3D) = [ panel.p1'; panel.p2'; panel.p3'; panel.p4'; panel.p1' ] 

collocation_point(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) / 4
panel_normal(panel :: Panel3D) = let p21 = panel.p2 .- panel.p1, 
                                     p41 = panel.p4 .- panel.p1, 
                                     p21_x_p41 = p21 × p41;
                                     p21_x_p41 / norm(p21_x_p41) end

weighted(x1, x2, wx) = (1 - wx) * x1 + wx * x2
weighted_point((x1, y1, z1), (x2, y2, z2), wx, wy, wz) = (weighted(x1, x2, wx), weighted(y1, y2, wy), weighted(z1, z2, wz))
quarter_point(p1, p2) = weighted_point(p1, p2, 1/4, 0, 1/4)
three_quarter_point(p1, p2) = weighted_point(p1, p2, 3/4, 0, 3/4)

bound_leg(p1, p2, p3, p4) = [ quarter_point(p1, p2); quarter_point(p4, p3) ]
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) .+ three_quarter_point(p4, p3) ) ./ 2

struct Line
    r1 :: SVector{3,Float64}
    r2 :: SVector{3,Float64}
end

function velocity(r, line :: Line, Γ, ε = 1e-6)
    r1 = r .- line.r1
    r2 = r .- line.r2
    r1_x_r2 = r1 × r2
    nr1, nr2 = norm(r1), norm(r2)

    any(<(ε), norm.([r1, r2, r1_x_r2])) ? [0,0,0] : Γ/(4π) * (1/nr1 + 1/nr2) * (r1_x_r2) / ( nr1 * nr2 + dot(r1, r2) )
    # Γ/(4π) * (r1_x_r2) / norm(r1_x_r2)^2 * dot(r1 - r2, r1 / nr1 - r2/nr2)
end

abstract type AbstractVortexArray <: Aircraft end

struct Horseshoe <: AbstractVortexArray
    vortex_lines :: Array{Line}
end

function Horseshoe(bound_leg :: Line, direction :: NTuple{3, Float64}, bound :: Float64 = 30.)
    r1, r2 = bound_leg.r1, bound_leg.r2
    vel = bound .* direction
    left_line = Line(r1 .+ vel, r1)
    right_line = Line(r2, r2 .+ vel)

    Horseshoe([ left_line, bound_leg, right_line ])
end

bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

vortex_collocation(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

horseshoe_lines(panel :: Panel3D, uniform :: Uniform3D) = Horseshoe(Line(bound_leg(panel)...), velocity(uniform) ./ uniform.mag)

struct VortexRing <: AbstractVortexArray
    vortex_lines :: Array{Line}
end

function vortex_ring(p1, p2, p3, p4)
    v1, v4 = SVector(quarter_point(p1, p2)...), SVector(quarter_point(p4, p3)...)
    v2, v3 = v1 .+ p2 .- p1, v4 .+ p3 .- p4
    [ v1, v2, v3, v4 ]
end 

vortex_ring(panel :: Panel3D) = vortex_ring(panel.p1, panel.p2, panel.p3, panel.p4)

function VortexRing(panel :: Panel3D)
    v1, v2, v3, v4 = vortex_ring(panel)

    bound_leg = Line(v1, v4)
    left_line = Line(v2, v1)
    right_line = Line(v3, v2)
    back_line = Line(v4, v3)

    VortexRing([ left_line, bound_leg, back_line, right_line ])
end

bound_leg_center(vortex_ring :: AbstractVortexArray) = let bound_leg = vortex_ring.vortex_lines[2]; (bound_leg.r1 .+ bound_leg.r2) / 2 end

bound_leg_vector(vortex_ring :: AbstractVortexArray) = let bound_leg = vortex_ring.vortex_lines[2]; (bound_leg.r2 .- bound_leg.r1) end

sum_vortices(r, vortex_lines :: Array{Line}, Γ) = sum(velocity(r, line, Γ) for line in vortex_lines)

velocity(r, vortex_ring :: AbstractVortexArray, Γ) = sum_vortices(r, vortex_ring.vortex_lines, Γ)

downwash_velocity(r, vortex_ring :: AbstractVortexArray, Γ) = let trailing_lines = [ first(vortex_ring.vortex_lines), last(vortex_ring.vortex_lines) ]; sum_vortices(r, trailing_lines, Γ) end

#---------------------------------Matrix setup--------------------------------------#

influence_coefficient(collocation_point, vortex_line :: AbstractVortexArray, panel_normal) = dot(velocity(collocation_point, vortex_line, 1.), panel_normal)

induced_coefficient(collocation_point, vortex_line :: AbstractVortexArray, panel_normal) = dot(downwash_velocity(collocation_point, vortex_line, 1.), panel_normal)

influence_matrix(colpoints, vortex_lines :: Array{<: AbstractVortexArray}, normals) = [ influence_coefficient(colpoint_i, horsie_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsie_j in vortex_lines ]

induced_matrix(colpoints, vortex_lines :: Array{<: AbstractVortexArray}, normals)  = [ induced_coefficient(colpoint_i, horsie_j, normal_i) for (colpoint_i, normal_i) in zip(colpoints, normals), horsie_j in vortex_lines ]

boundary_condition(normals, velocity) = - [ dot(velocity, normal) for normal in normals ]

function accumap(f, n, xs)
    data = [ xs ]
    for i = 1:n
        ys = map(f, xs)
        data = [ data..., ys ]
        xs = ys
    end
    return hcat(data...)
end

function wake_panel(panel :: Panel3D, uniform :: Uniform3D, bound = 1.)
    wake = bound .* velocity(uniform) ./ uniform.mag
    Panel3D(panel.p2, panel.p2 .+ wake, panel.p3 .+ wake, panel.p3)
end

make_wake(last_panels :: Array{Panel3D}, uniform :: Uniform3D, wake_length, wake_num) = (permutedims ∘ accumap)(panel -> wake_panel(panel, uniform, wake_length / wake_num), wake_num, last_panels)


#-------------------------Force evaluations------------------------------------#

# Trefftz plane evaluations
# downwash(xi, xj, zi, zj) = -1/(2π) * (xj - xi) / ( (zj - zi)^2 + (xj - xi)^2 )
# trefftz_plane(horseshoes :: Array{Horseshoes}) = [ [ horseshoe.vortex_lines[1].r1 for horseshoe in horseshoes ]; 
#                                                    [ horseshoe.vortex_lines[3].r2 for horseshoe in horseshoes ] ]
# function downwash(horseshoes :: Array{Horseshoes}) 
#     coords = trefftz_plane(horseshoes)
#     carts = product(coords, coords)
#     [ [ i == j ? 0 : downwash(x1, x2, z1, z2) for (i, (x1,y1,z1)) in enumerate(coords) ] for (j, (x2,y2,z2)) in enumerate(coords) ] 
# end

# Horseshoe forces and moments
function kutta_force(Γs, vortices :: Array{<: AbstractVortexArray}, uniform :: Uniform3D, ρ = 1.225)
    Γ_rings = zip(Γs, vortices)
    [ sum(velocity(uniform) .+ velocity(bound_leg_center(vortex_focus), vortex, Γ) for (Γ, vortex) in Γ_rings) × bound_leg_vector(vortex_focus) * ρ * Γ_focus for (Γ_focus, vortex_focus) in Γ_rings ] 
end

# forces(vortex_legs :: Array{<: AbstractVortexArray}, Γs, uniform :: Uniform3D, ρ) = [ kutta_force(Γ, vortex_ring, uniform, ρ) for (Γ, vortex_ring) in zip(Γs, vortex_legs) ]

moments(vortex_legs :: Array{<: AbstractVortexArray}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) in zip(forces, vortex_legs) ]

# Aerodynamic coefficients
force_coefficient(F, ρ, V, S) = F / (0.5 * ρ * V^2 * S)
moment_coefficient(M, ρ, V, S, ref_length) = M / (0.5 * ρ * V^2 * S * ref_length)

#------------------------Case setup and solution--------------------------#

function solve_horseshoes(panels :: Array{Panel3D}, uniform :: Uniform3D, r_ref = (0.25, 0, 0), ρ = 1.225) 
    
    vel = velocity(uniform)
    horseshoes = [ horseshoe_lines(panel, uniform) for panel in panels ][:]
    colpoints = vortex_collocation.(panels)[:]
    normals = panel_normal.(panels)[:]
    horsies = bound_leg.(panels)[:]

    inf = influence_matrix(colpoints, horseshoes, normals)
    
    boco = boundary_condition(normals, vel)

    Γs = inf \ boco

    println(Γs)

    # w_inds = induced_matrix(colpoints, horseshoes, normals) * Γs
    # Δys = (norm ∘ bound_leg_vector).(horseshoes)
    # lols = zip(Γs, Δys, w_inds)
    # lift = [ ρ*uniform.mag*Γj*Δyj for (Γj, Δyj, w_indj) in lols ]
    # drag = [ -ρ*uniform.mag*w_indj*Γj*Δyj for (Γj, Δyj, w_indj) in lols ]

    # # Force computations
    geom_forces = kutta_force(Γs, horseshoes, uniform, ρ)
    geom_moments = moments(horseshoes, geom_forces, r_ref)

    force, moment = sum(geom_forces), sum(geom_moments)
    drag = dot(force, vel ./ uniform.mag)
    
    stable_forces = stability_axes(force, uniform)
    stable_moments = stability_axes([ moment[1], -moment[2], moment[3] ], uniform)

    geom_forces, geom_moments, stable_forces, stable_moments, drag, Γs, horseshoes
    # lift, drag
end

function streamlines(point, uniform :: Uniform3D, horseshoes, Γs, length, num_steps = 100)
    streamlines = [point]
    for i in 1:num_steps
        update = velocity(uniform) .+ sum(velocity(streamlines[end], horseshoe, Γ) for (horseshoe, Γ) in zip(horseshoes, Γs))
        streamlines = [ streamlines..., streamlines[end] .+ (update/norm(update) * length / num_steps)  ]
    end
    streamlines
end

function solve_vortex_rings(panels :: Array{Panel3D}, uniform :: Uniform3D, wake_num = 10, wake_length = 30., r_ref = (0.25, 0, 0)) 
    
    vel = velocity(uniform)

    first_panels = panels[1,:]
    last_panels = panels[end,:]

    wake_panels = make_wake(last_panels, uniform, wake_length, wake_num)
    woke_panels = [ panels;
                    wake_panels ]

    vortex_legs_2D = VortexRing.(woke_panels)
    vortex_legs = vortex_legs_2D[:]
    colpoints = vortex_collocation.(woke_panels)[:]
    normals = panel_normal.(woke_panels)[:]
    horsies = bound_leg.(woke_panels)[:]

    inf = influence_matrix(colpoints, vortex_legs, normals)
    boco = boundary_condition(normals, vel)

    Γs =  inf \ boco

    # # Force computations
    Γs = reshape(Γs, (length(woke_panels[:,1]), length(woke_panels[1,:])))[1:end-wake_num,:]
    vortex_legs_2D = vortex_legs_2D[1:end-wake_num,:]

    le_forces = kutta_force(Γs[1,:], vortex_legs_2D[1,:], uniform, ρ)
    trailing_forces = kutta_force(map(-, Γs[3:end,:], Γs[2:end-1,:]), vortex_legs_2D[2:end-1,:], uniform, ρ)

    geom_forces = [ permutedims(le_forces); trailing_forces ][:]
    
    geom_moments = moments(vortex_legs, geom_forces, r_ref)

    pressures = geom_forces ./ area.(panels[:])

    force, moment = sum(geom_forces), sum(geom_moments)
    drag = dot(force, vel ./ uniform.mag)
    
    stable_forces = stability_axes(force, uniform)
    stable_moments = stability_axes([ -moment[1], moment[2], -moment[3] ], uniform)

    geom_forces, geom_moments, stable_forces, stable_moments, pressures, drag, woke_panels[:]
end

function print_dynamics(force, moment, drag, V, S, b, c, ρ = 1.225)
    CDi = force_coefficient(drag, ρ, V, S)
    CY = force_coefficient(force[2], ρ, V, S)
    CL = force_coefficient(force[3], ρ, V, S)

    Cl = moment_coefficient(moment[1], ρ, V, S, b)
    Cm = moment_coefficient(moment[2], ρ, V, S, c)
    Cn = moment_coefficient(moment[3], ρ, V, S, b)

    println("Total Force: $force N")
    println("Total Moment: $moment N-m")
    println("Lift Coefficient (CL): $CL")
    println("Drag Coefficient (CDi): $CDi")
    println("Side Force Coefficient (CY): $CY")
    println("Lift-to-Drag Ratio (L/D): $(CL/CDi)")
    println("Rolling Moment Coefficient (Cl): $Cl")
    println("Pitching Moment Coefficient (Cm): $Cm")
    println("Yawing Moment Coefficient (Cn): $Cn")
end

end