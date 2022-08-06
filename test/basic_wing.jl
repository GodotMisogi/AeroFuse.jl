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
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

## Freestream velocity
alpha = 5.0
beta = 0.0
fs = Freestream(alpha, beta, zeros(3))
V∞ = 15.

V = V∞ * velocity(fs)

## Wake panels
wake_pans = wake_panel.(eachcol(surf_pans), 100., alpha, beta)

## Aerodynamic Influence Coefficient matrix
# Foil panel interaction
npanscd, npanssp = size(surf_pans)
npans = prod([npanscd, npanssp])
AIC_ff = doublet_matrix(surf_pans[:], surf_pans[:])

# Wake panel interaction
AIC_wf = doublet_matrix(surf_pans[:], wake_pans[:])

# Kutta condition
AIC = [	AIC_ff				AIC_wf		;
							I(npanssp)	]

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