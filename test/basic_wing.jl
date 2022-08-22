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
wing = WingSection(
    area      = 19.0, # Projected area
    aspect    = 8.3, # Aspect ratio
    dihedral  = 0.0, # Dihedral angle (deg)
    sweep     = 0.0, # Sweep angle (deg)
    w_sweep   = 0.25, # Sweep location normalized to chords ∈ [0,1]
    taper     = 1.0, # Taper ratio
    root_foil = naca4(0,0,1,2), # Root airfoil
    tip_foil  = naca4(0,0,1,2), # Tip airfoil
    position  = [0.0, 0, 0.0], # Location (m)
    angle     = 0., # Angle of incidence (deg)
    axis      = [0., 1., 0.0], # Spanwise direction
    symmetry  = true # Symmetric about x-z plane
)

# Meshing
wing_mesh = WingMesh(wing, [40], 40; span_spacing = Uniform())

surf_pts  = surface_coordinates(wing_mesh)
left_edge = (surf_pts[:,1] + reverse(surf_pts[:,1])) / 2
right_edge = (surf_pts[:,end] + reverse(surf_pts[:,end])) / 2
surf_pts = [left_edge surf_pts right_edge]
surf_pans = make_panels(surf_pts)

# Freestream velocity
α = 6.0; β = 0.0; Umag = 15.
fs = Freestream( α, β, zeros(3))

##
@time prob = solve_system_neumann(surf_pans, Umag, fs, 5)
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

##
using Plots
plotlyjs()
Plots.plot()
plt = plot_panels(prob.surface_panels)
[f = Plots.plot!(p, color = :blue, zlim=(-2.5,2.5)) for p in plt]
plt = plot_panels(prob.wake_panels)
[f = Plots.plot!(p, color = :grey, zlim=(-2.5,2.5)) for p in plt]
Plots.plot!()