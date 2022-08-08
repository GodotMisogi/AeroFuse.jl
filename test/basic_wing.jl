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
                   taper      = 1.0,
                   span       = 4.0,
                   dihedral   = 0.0,
                   sweep      = 0.0)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 0.0
β = 0.0
Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
system = solve_system(surf_pans, Umag, fs, 1.0e5)

##
vxs, vys = surface_velocities(system)
@time cls, cps = surface_coefficients(system);
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots
plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-2.0, 3.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

##
