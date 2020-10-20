## 
using Revise
includet("../src/AeroMDAO.jl")
includet("../src/FoilParametrization.jl")
includet("../src/MathTools.jl")
includet("../src/Geometry.jl")

##
using .AeroMDAO: Foil, HalfWing, Wing, projected_area, span, mean_aerodynamic_chord, horseshoe_points, horseshoe_collocation, wing_coords, wing_sections, make_panels
using .FoilParametrization: read_foil
using .MathTools: linspace
using .Geometry: Point2D, Point3D
using DelimitedFiles
using Rotations

## Wing section setup
foilpath = "airfoil_database/ys930.dat"

num_secs = 6

coords = read_foil(foilpath)
foils = [ coords for n in 1:num_secs ]
airfoils = Foil.(foils)

wing_chords = [2, 2, 1.5, 1, 0.8, 0.7, 0.5]
wing_twists = fill(1, num_secs + 1)
wing_spans = fill(1, num_secs)
wing_dihedrals = 0.5 * [0, 5, 10, 10, 15, 20]
wing_sweeps = [0, 5, 10, 10, 15, 20]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
println("Span: ", span(wing), " m")
println("Area: ", projected_area(wing), " m²")
println("MAC: ", mean_aerodynamic_chord(wing), " m")
wing_lead, wing_trail = wing_coords(wing_right)
println("Leading")
writedlm(stdout, wing_lead)
println("Trailing")
writedlm(stdout, wing_trail)
wing_secs = wing_sections(wing)
panels = make_panels(wing)

## Horizontal tail section setup
tail_secs = 1

tail_chords = [0.8, 0.4]
tail_twists = fill(0, tail_secs + 1)
tail_spans = fill(1, tail_secs)
tail_dihedrals = [5]
tail_sweeps = [30]

tail_right = HalfWing(airfoils, tail_chords, tail_spans, tail_dihedrals, tail_sweeps, tail_twists)
tail = Wing(tail_right, tail_right)
println("Span: ", span(tail), " m")
println("Area: ", projected_area(tail), " m²")
println("MAC: ", mean_aerodynamic_chord(tail), " m")
htaily = wing_coords(tail) .+ [5 0 0]
writedlm(stdout, htaily)

## Vertical tail setup
vtail_secs = 1

vtail_chords = [0.8, 0.4]
vtail_twists = fill(0, vtail_secs + 1)
vtail_spans = fill(0.4, vtail_secs)
vtail_dihedrals = [0.0]
vtail_sweeps = [60]

vtail_right = HalfWing(airfoils, vtail_chords, vtail_spans, vtail_dihedrals, vtail_sweeps, vtail_twists)
println("Span: ", span(vtail_right), " m")
println("Area: ", projected_area(vtail_right), " m²")
println("MAC: ", mean_aerodynamic_chord(vtail_right), " m")
vtaily = wing_coords(vtail_right) * AngleAxis{Float64}(π/2, 1, 0, 0)' .+ [5 0 0]

## Plotting
using Plots, LaTeXStrings
plotlyjs()

plot(wingy[:,1], wingy[:,2], wingy[:,3], label = "Wing")
plot!(htaily[:,1], htaily[:,2], htaily[:,3], label = "Horizontal Tail")
plot!(vtaily[:,1], vtaily[:,2], vtaily[:,3], label = "Vertical Tail",
      xaxis = L"x", yaxis = L"y", zaxis = L"z", zlim = (-0.1, 5), aspect_ratio=:equal)

plot!(wing_secs[:,1], wing_secs[:,2], wing_secs[:,3], color = "black", label = nothing)

## Horseshoe testing
pts = horseshoe_points(wing_right)
collocs = horseshoe_collocation(wing_right)
println("Horseshoe Points: ", pts)
println("Horseshoe Collocations: ", collocs)