## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/MathTools.jl")

##
using .AeroMDAO
using .MathTools: linspace, tuparray
using DelimitedFiles
using Rotations

## Wing section setup
foilpath = "airfoil_database/ys930.dat"

num_secs = 3

coords = read_foil(foilpath)
foils = [ coords for i in 1:num_secs ]
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
test_panels = mesh_wing(wing, spanwise_panels = 5, chordwise_panels = 5)
cam_panels = mesh_cambers(wing_right, spanwise_panels = 1, chordwise_panels = 2)
# wing_panels = make_panels(wing_right, spanwise_panels = 5, chordwise_panels = 1)

## Panel case
ρ = 1.225
uniform = Uniform(10.0, 5.0, 0.0)
@time lift, drag = solve_case(wing_panels, uniform, ρ)

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
# wing_pts = horseshoe_vortex.(wing_panels)
# wing_collocs = horseshoe_collocation.(cam_panels)
pan_coords = (tuparray ∘ panel_coords).(test_panels)
cam_coords = (tuparray ∘ panel_coords).(cam_panels)
# spans = [ pt[2] for pt in wing_collocs ];

##
plot(spans, cl)
plot(spans, cdi)

##
plot(xaxis = "x", yaxis = "y", zaxis = "z", aspect_ratio = :equal, zlim = (-0.5, 5.0))
plot!.(pan_coords, color = :black,label = :none)
plot!.(cam_coords, color = :grey,label = :none)

scatter!(wing_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")

# plot!.(wing_pts, c = :black, label = :none)

# wing_lines = [ horseshoe_lines(panel, uniform) for panel in wing_panels ]

# lines = [ [ tuparray([line.r1'; line.r2']) for line in horseshoe.vortex_lines ] for horseshoe in wing_lines ]
# [ [ plot!(line, c = :darkblue, label = :none) for line in horseshoe ] for horseshoe in lines ]

gui();