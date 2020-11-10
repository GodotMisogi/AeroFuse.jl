## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")

##
using StaticArrays
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

wing_chords = [1, 0.5, 0.1]
wing_twists = [4, 2, -2]
wing_spans = [2, 0.2]
wing_dihedrals = [0, 30]
wing_sweeps = [15, 45]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
ρ = 1.225
uniform = Uniform(10.0, 5.0, -3.0)
@time horseshoe_panels, horseshoes, Γs = solve_case(wing, uniform, span_num = 15, chord_num = 5);

## Panel method: TO DO
camber_panels = mesh_cambers(wing, span_num = 5, chord_num = 10)
wing_panels = mesh_wing(wing, span_num = 5, chord_num = 10)

##
using Plots, LaTeXStrings
##
plotlyjs()

##
wing_coords = plot_panels(wing_panels)[:]
horseshoe_coords = plot_panels(horseshoe_panels)[:]
streams = plot_streamlines.(streamlines(uniform, horseshoe_panels, horseshoes, Γs, 5, 100))

##
plot(xaxis = "x", yaxis = "y", zaxis = "z", aspect_ratio = :equal, zlim = (-0.5, 5.0), size=(1280, 720))
plot!.(wing_coords, color = :black, label = :none)
plot!.(horseshoe_coords, color = :grey, label = :none)
plot!.(streams, color = :green, label = :none)

gui();