##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   tip_foil   = naca4(0,0,1,2),
                   root_chord = 1.0,
                   taper      = 0.6,
                   span       = 5.0,
                   dihedral   = 0.0,
                   sweep      = 0.0)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)
surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

# Freestream velocity
α = 10.0; β = 0.0; Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
prob = solve_system_neumann(surf_pans, Umag, fs, 50)
φs = permutedims(reshape(prob.singularities[1:360], size(prob.surface_panels,2), size(prob.surface_panels,1)))

##
vs = surface_velocities(prob)
@time cls, cps = surface_coefficients(prob, wing);
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots
plotly()
plt_surfs = plot_panels(surf_pans)

p = Plots.plot(; xlabel="x",ylabel="y",zlabel="z", aspect_ratio=:equal, grid=:true, zlim=(-2.0, 3.0))

for pan in plt_surfs
    p = Plots.plot!(pan, color = :grey, aspect_ratio=:equal)
end

plt_surfs = plot_panels(prob.wake_panels)
for pan in plt_surfs
    p = Plots.plot!(pan, color = :blue, aspect_ratio=:equal)
end

plot!(p)
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

##
