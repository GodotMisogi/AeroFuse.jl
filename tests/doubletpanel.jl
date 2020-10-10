include("../src/PanelMethods.jl")
include("../src/FoilParametrization.jl")

using .PanelMethods: DoubletPanel2D, Uniform2D, grid_data, aero_coefficients
using .FoilParametrization: read_foil, cosine_foil
using BenchmarkTools
using Plots
plotly()

airfoil = cosine_foil(read_foil("airfoil_database/s1223.dat"), n = 100)

panels = reverse([ DoubletPanel2D((xs, ys), (xp, yp)) for (xs, ys, xp, yp) ∈ (collect ∘ eachrow)([airfoil[2:end,:] airfoil[1:end-1,:]]) ], dims=1)

cls = []
for alpha in 0:15
    uniform = Uniform2D(5.0, alpha)

    @time cl, kcl, cps = aero_coefficients(panels, uniform)
    println("Angle: $alpha, Lift Coefficient: $cl, $kcl")

    append!(cls, cl)
end
plot(0:15, cls)