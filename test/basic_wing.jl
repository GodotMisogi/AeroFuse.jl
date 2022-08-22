##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO
using LinearAlgebra
using GLMakie
Makie.inline!(false)
# using Makie.GeometryBasics

## Geometry
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   tip_foil   = naca4(0,0,1,2),
                   root_chord = 1.0,
                   taper      = 0.4,
                   span       = 5.0,
                   dihedral   = 5.0,
                   sweep      = 25.0)

# Meshing
wing_mesh = WingMesh(wing, [20], 20; span_spacing = Uniform())

surf_pts  = surface_coordinates(wing_mesh)
left_edge = (surf_pts[:,1] + reverse(surf_pts[:,1])) / 2
right_edge = (surf_pts[:,end] + reverse(surf_pts[:,end])) / 2
surf_pts = [left_edge surf_pts right_edge]
surf_pans = make_panels(surf_pts)

# Freestream velocity
α = 8.0; β = 0.0; Umag = 15.
fs = Freestream( α, β, zeros(3))

##
@time prob = solve_system_neumann(surf_pans, Umag, fs, 50)
@time vs = surface_velocities(prob);
@time CL, CD, CP = surface_coefficients(prob, wing);
println("Σᵢ CLᵢ: $CL")
println("Σᵢ CDᵢ: $CD")

## Plotting
f = Figure()
ax = LScene(f[1, 1])
pts, cnt, C = plot_scalar_field(surf_pans, CP)
poly!(pts, cnt; color=C[:], strokewidth=0.1)
current_figure()

