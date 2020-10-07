include("../src/AeroMDAO.jl")
include("../src/FoilParametrization.jl")

using .AeroMDAO: Point3D, Point2D, Foil, HalfWing, Wing, projected_area, span, mean_aerodynamic_chord
using .FoilParametrization: read_foil, linspace
using PyPlot

foilpath = "airfoil_database/ys930.dat"

# Wing section setup
num_secs = 5
xs = zeros(num_secs)
ys = linspace(0, 2, num_secs)
zs = zeros(num_secs)

coords = read_foil(foilpath)
foils = [ coords for n in 1:num_secs ]
airfoils = Foil.(foils)

chords = [2, 1.5, 1, 0.5, 0.2, 0.1]
twists = zeros(num_secs + 1)
spans = repeat([1.0], num_secs)
dihedrals = repeat([5.0], num_secs)
sweeps = repeat([30.0], num_secs)

left = HalfWing(airfoils, chords, spans, dihedrals, sweeps, twists)
wing = Wing(left, left)
println("Span: ", span(wing), " m")
println("Area: ", projected_area(wing), " m²")
println("MAC: ", mean_aerodynamic_chord(wing), " m")

# wing_loc = Point3D{Float64}(0.0, 0.0, 0.0)
# wing = Wing(wing_loc, secs)

# println(projected_area(wing))