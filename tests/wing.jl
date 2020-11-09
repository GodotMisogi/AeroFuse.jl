## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: read_foil, foil_camthick, camthick_foil, cosine_foil, kulfan_CST, naca4
using .AeroMDAO
using .MathTools: linspace, tuparray, tupvector
using DelimitedFiles
using Rotations

## Wing section setup
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
alphas = [alpha_u alpha_l]
dzs = (1e-4, 1e-4)
cst_foil = kulfan_CST(alphas, dzs, 0.2)

num_secs = 3
foils = [ cst_foil for i in 1:num_secs ]
airfoils = Foil.(foils)

wing_chords = [0.18, 0.16, 0.08]
wing_twists = [2, 0, -2]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0, 5]
wing_sweeps = [0, 30]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing_right)

## Assembly
ρ = 1.225
uniform = Uniform(10.0, 5.0, 0.0)
V = uniform.mag
S = projected_area(wing_right)
b = span(wing_right)
c = mean_aerodynamic_chord(wing_right);

## Horseshoe method
horseshoe_panels = mesh_horseshoes(wing, spanwise_panels = 5, chordwise_panels = 5)
@time forces, moments, stable_force, stable_moment, drag, gammas, horseshoes = solve_horseshoes(horseshoe_panels, uniform)

# print(forces)
# Post-processing
print_dynamics(stable_force, stable_moment, drag, V, S, b, c, ρ)

## Vortex lattice method
# camber_panels = mesh_cambers(wing_right, spanwise_panels = 5, chordwise_panels = 5)
# @time forces, moments, stable_force, stable_moment, pressures, drag, camber_panels = solve_vortex_rings(camber_panels, uniform, 10, 30, (0.25, 0, 0))

## Post-processing
# print_dynamics(stable_force, stable_moment, drag, V, S, b, c, ρ)


## Panel method: TO DO
wing_panels = mesh_wing(wing, spanwise_panels = 5, chordwise_panels = 10)

##
using Plots, LaTeXStrings
##
plotlyjs()


## Wing
horseshoe_collocs = vortex_collocation.(horseshoe_panels)[:]
# camber_collocs = vortex_collocation.(camber_panels)[:]

## Streamlines
streams = [ tupvector(streamlines(SVector(hs), uniform, horseshoes, gammas, 5, 100)) for hs in horseshoe_collocs ]

##
spans = [ pt[2] for pt in horseshoe_collocs ];
plot(spans, lift, label = L"C_L")
plot!(spans, drag, label = L"C_{D_i}")

##
# wing_coords = (tuparray ∘ panel_coords).(wing_panels)[:]
horseshoe_coords = tuparray.(panel_coords.(horseshoe_panels)[:])
# camber_coords = (tuparray ∘ panel_coords).(camber_panels)[:]
# vortex_rings = (tupvector ∘ vortex_ring).(camber_panels)

##
plot(xaxis = "x", yaxis = "y", zaxis = "z", aspect_ratio = :equal, zlim = (-0.5, 5.0), size=(1280, 720))
plot!.(wing_coords, color = :black, label = :none)
plot!.(horseshoe_coords, color = :grey, label = :none)
plot!.(streams, color = :lightblue, label = :none)
# plot!.(camber_coords, color = :grey, label = :none)
# plot!.(vortex_rings, color = :black, label = :none)

# gui();

scatter!(horseshoe_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")
# scatter!(camber_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")

##
wing_lines = [ horseshoe_lines(panel, uniform) for panel in horseshoe_panels ]

lines = [ [ tuparray([line.r1'; line.r2']) for line in horseshoe.vortex_lines ] for horseshoe in wing_lines ]
[ [ plot!(line, c = :darkblue, label = :none) for line in horseshoe ] for horseshoe in lines ]

gui();