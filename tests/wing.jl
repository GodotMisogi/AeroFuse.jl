## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/FoilParametrization.jl")
includet("../src/MathTools.jl")
# includet("../src/LiftingLine.jl")

##
using .AeroMDAO
# using .LiftingLine
using .FoilParametrization: read_foil
using .MathTools: linspace
using DelimitedFiles
using Rotations

## Wing section setup
foilpath = "airfoil_database/ys930.dat"

num_secs = 6

coords = read_foil(foilpath)
foils = [ coords for n in 1:num_secs ]
airfoils = Foil.(foils)

wing_chords = [2, 2, 1.5, 1, 0.8, 0.7, 0.5]
wing_twists = fill(0, num_secs + 1)
wing_spans = fill(1, num_secs)
wing_dihedrals = 0.0 * [0, 5, 10, 10, 15, 20]
wing_sweeps = [0, 5, 10, 10, 15, 20]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing)

## Horizontal tail section setup
tail_secs = 1

htail_chords = [0.8, 0.4]
htail_twists = fill(0, tail_secs + 1)
htail_spans = fill(1, tail_secs)
htail_dihedrals = [0]
htail_sweeps = [30]
htail_location = [5 0 0]

htail_right = HalfWing(airfoils, htail_chords, htail_spans, htail_dihedrals, htail_sweeps, htail_twists)
htail = Wing(htail_right, htail_right)
print_info(htail)


## Vertical tail setup
vtail_secs = 1

vtail_chords = [0.4, 0.2]
vtail_twists = fill(0, vtail_secs + 1)
vtail_spans = fill(0.3, vtail_secs)
vtail_dihedrals = [0.0]
vtail_sweeps = [60]


vtail = HalfWing(airfoils, vtail_chords, vtail_spans, vtail_dihedrals, vtail_sweeps, vtail_twists)
print_info(vtail)

vtail1_angle = AngleAxis{Float64}(π/4, 1, 0, 0)
vtail1_location = [5 + sin(π/6) 1 0]

vtail2_angle = AngleAxis{Float64}(3π/4, 1, 0, 0)
vtail2_location = [5 + sin(π/6) -1 0]

## Panelling
wing_panels = make_panels(wing)
htail_panels = make_panels(htail, translation = htail_location)
vtail1_panels = make_panels(vtail, rotation = vtail1_angle, translation = vtail1_location)
vtail2_panels = make_panels(vtail, rotation = vtail2_angle, translation = vtail2_location)

## Horseshoe testing
wing_pts = horseshoe_vortex.(wing_panels)
wing_collocs = horseshoe_collocation.(wing_panels)

htail_pts = horseshoe_vortex.(htail_panels)
htail_collocs = horseshoe_collocation.(htail_panels)

vtail1_pts = horseshoe_vortex.(vtail1_panels)
vtail1_collocs = horseshoe_collocation.(vtail1_panels)

vtail2_pts = horseshoe_vortex.(vtail2_panels)
vtail2_collocs = horseshoe_collocation.(vtail2_panels)

wing_lines = horseshoe_lines.(wing_panels)

## Solution for wing
uniform = Uniform3D(5.0, 4.0, 0.0)
@time lift, drag = solve_case(wing_panels, uniform)

V = uniform.mag
ρ = 1.225
S = projected_area(wing)
cl = lift/(ρ/2 * V^2 * S)
cdi = drag/(ρ/2 * V^2 * S)
println("Lift: $lift, Drag: $drag")
println("Lift Coefficient: $cl, Drag Coefficient: $cdi")

## Plotting
using Plots, LaTeXStrings
plotlyjs()

wing_plot = plot_setup(coordinates(wing))
htail_plot = plot_setup(coordinates(htail, translation = htail_location))
vtail1_plot = plot_setup(coordinates(vtail, rotation = vtail1_angle, translation = vtail1_location))
vtail2_plot = plot_setup(coordinates(vtail, rotation = vtail2_angle, translation = vtail2_location))

wing_sects = plot_setup.(sections(wing))
htail_sects = plot_setup.(sections(htail, translation = htail_location))
vtail1_sects = plot_setup.(sections(vtail, rotation = vtail1_angle, translation = vtail1_location))
vtail2_sects = plot_setup.(sections(vtail, rotation = vtail2_angle, translation = vtail2_location))

# Objects
plot(wing_plot, label = "Wing", fill = :blue)
plot!(htail_plot, label = "Horizontal Tail")
plot!(vtail1_plot, label = "Vertical Tail 1")
plot!(vtail2_plot, label = "Vertical Tail 2",
      xaxis = L"x", yaxis = L"y", zaxis = L"z", zlim = (-0.1, 5), aspect_ratio=:equal)

# Sections
plot!.(wing_sects, color = "black", label = nothing)
plot!.(htail_sects, color = "orange", label = nothing)
plot!.(vtail1_sects, color = "brown", label = nothing)
plot!.(vtail2_sects, color = "brown", label = nothing)

# Points
plot!.(wing_pts, c = :black, label = nothing)
scatter!(wing_collocs, c = :grey, markersize = 1, label = "Wing Collocation Points")

plot!.(htail_pts, c = :black, label = nothing)
scatter!(htail_collocs, c = :grey, markersize = 1, label = "Horizontal Tail Collocation Points")

plot!.(vtail1_pts, c = :black, label = nothing)
scatter!(vtail1_collocs, c = :grey, markersize = 1, label = "Vertical Tail 1 Collocation Points")

plot!.(vtail2_pts, c = :black, label = nothing)
scatter!(vtail2_collocs, c = :grey, markersize = 1, label = "Vertical Tail 2 Collocation Points")

gui();
