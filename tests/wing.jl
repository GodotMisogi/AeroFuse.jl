## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: read_foil, foil_camthick, camthick_foil, cosine_foil, kulfan_CST, naca4
using .AeroMDAO
using .MathTools: linspace, tuparray
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

wing_chords = [1, 0.5, 0.2]
wing_twists = [5, 2, 0]
wing_spans = [2, 0.5]
wing_dihedrals = [0, 45]
wing_sweeps = [15, 30]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing_right)

## Assembly
wing_panels = mesh_wing(wing, spanwise_panels = 5, chordwise_panels = 20)
camber_panels = mesh_cambers(wing, spanwise_panels = 1, chordwise_panels = 20)
horseshoe_panels = mesh_horseshoes(wing_right, spanwise_panels = 5, chordwise_panels = 20)

## Panel case
ρ = 1.225
uniform = Uniform(10.0, 5.0, 0.0)
@time lift, drag = solve_case(horseshoe_panels, uniform, ρ)

V = uniform.mag
S = projected_area(wing_right)
cl = lift/(ρ/2 * V^2 * S)
cdi = drag/(ρ/2 * V^2 * S)
println("Lift: $(sum(lift)) N")
println("Drag: $(sum(drag)) N")
println("Lift Coefficient: $(sum(cl)), Drag Coefficient: $(sum(cdi))")
println("Lift-to-Drag Ratio (L/D): $(sum(cl)/sum(cdi))")


##
using Plots
##
plotlyjs()


## Wing
wing_collocs = horseshoe_collocation.(horseshoe_panels)

##
# spans = [ pt[2] for pt in wing_collocs ];
# plot(spans, cl)
# plot(spans, cdi)

##
pan_coords = (tuparray ∘ panel_coords).(wing_panels)
wing_coords = (tuparray ∘ panel_coords).(horseshoe_panels)
cam_coords = (tuparray ∘ panel_coords).(camber_panels)


##
plot(xaxis = "x", yaxis = "y", zaxis = "z", aspect_ratio = :equal, zlim = (-0.5, 5.0), size=(800, 600))
plot!.(pan_coords, color = :black,label = :none)
plot!.(cam_coords, color = :grey,label = :none)
plot!.(wing_coords, color = :grey,label = :none)

scatter!(wing_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")

##
wing_lines = [ horseshoe_lines(panel, uniform) for panel in horseshoe_panels ]

lines = [ [ tuparray([line.r1'; line.r2']) for line in horseshoe.vortex_lines ] for horseshoe in wing_lines ]
[ [ plot!(line, c = :darkblue, label = :none) for line in horseshoe ] for horseshoe in lines ]

gui();