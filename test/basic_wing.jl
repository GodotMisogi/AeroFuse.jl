##
using Revise
using AeroMDAO
using LinearAlgebra

## Geometry
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   tip_foil   = naca4(0,0,1,2),
                   root_chord = 1.0,
                   taper      = 0.6,
                   span       = 4.0,
                   dihedral   = 0.0,
                   sweep      = 0.0)

## Meshing
wing_mesh = WingMesh(wing, [15], 30)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)
wake_pans = wake_panels(surf_pans)

## Freestream velocity
fs = Freestream(alpha = 5.0)
V∞ = 15.

V = V∞ * velocity(fs)

## Aerodynamic Influence Coefficient matrix
# Foil panel interaction
npanscd, npanssp = size(surf_pans)
npans = prod([npanscd, npanssp])
AIC_ff = doublet_matrix(surf_pans[:], surf_pans[:])

# Wake panel interaction


# Boundary condition (no sources yet)
RHS = [ dot(V∞, pt) for p_j in ??? ]

# Solve
φ = AIC \ RHS

## Plotting
using Plots

plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.2, 1.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)