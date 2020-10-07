include("../src/PanelMethods.jl")
include("../src/FoilParametrization.jl")

using .PanelMethods: DoubletPanel2D, Uniform2D, grid_data, aero_coefficients
using .FoilParametrization: read_foil, cosine_foil
using BenchmarkTools
using PyPlot

coords = cosine_foil(read_foil("airfoil_database/s1223.dat"), n = 5)

panels = reverse([ DoubletPanel2D((xs, ys), (xp, yp)) for (xs, ys, xp, yp) ∈ (collect ∘ eachrow)([coords[1:end-1,:] coords[2:end,:]]) ], dims=1)

uniform = Uniform2D(5.0, 5.0)

@time (cl, kcl, cps) = aero_coefficients(panels, uniform)
println("Lift Coefficient: $cl, $kcl")