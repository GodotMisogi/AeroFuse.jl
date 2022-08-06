##
using Pkg
Pkg.activate(".")
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

foil_pans = @view permutedims(surf_pans)[:]

## Freestream velocity
alpha = 5.0
beta = 0.0
Umag = 15
fs = Freestream(alpha, beta, zeros(3))
V∞ = Umag * velocity(fs)

## Wake panels
wake_pans = wake_panel.(eachcol(surf_pans), 1.0e5, alpha, beta)

## Aerodynamic Influence Coefficient matrix
npancd, npansp = size(surf_pans)
npanf = npancd * npansp
npanw = npansp

# Foil panel interaction
AIC_ff = doublet_matrix(foil_pans, foil_pans)

# Wake panel interaction
AIC_wf = doublet_matrix(foil_pans, wake_pans)

# Kutta condition
AIC_ff[:,1:npanw] 			+= AIC_wf
AIC_ff[:,end-npanw+1:end] 	-= AIC_wf

## Boundary condition (no sources yet)
boco = [ dot(V∞, pt) for pt in collocation_point.(foil_pans) ]

## Solve
φ = AIC_ff \ boco

##
problem = DoubletSourceSystem3D(AIC_ff, boco, φ, foil_pans, wake_pans, fs, Umag)

## Plotting
using Plots

plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.2, 1.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)