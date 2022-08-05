##
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

## Freestream velocity
fs = Freestream(alpha = 5.0)
V∞ = 15.

V = V∞ * velocity(fs)

## Template:

# Aerodynamic Influence Coefficient matrix
npan = len(surf_pans)
AIC_ff = doublet_matrix(panels, panels)

for i=1:npan
	tr = get_transformation(surf_pans[i])
	for j=npan
		point = collocation_point(surf_pans[j])
		AIC_ff[j, i] = ifelse(
			surf_pans[i] == surf_pans[j],
			0.5,
			quadrilateral_doublet_potential(1., panel, point)
		)
end

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