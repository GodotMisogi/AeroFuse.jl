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
                   sweep      = 20.0)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 0
β = 0.0
Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
prob = solve_system(surf_pans, Umag, fs, 1)














##
vxs, vys = surface_velocities(prob)
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
