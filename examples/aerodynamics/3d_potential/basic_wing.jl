##
using AeroFuse
using LinearAlgebra

## Geometry
wing = WingSection(
    root_foil  = naca4(0,0,1,2),
    tip_foil   = naca4(0,0,1,2),
    taper      = 0.6,
    aspect     = 6.0,
    dihedral   = 0.0,
    sweep      = 0.0,
    symmetry   = true
)

## Meshing
wing_mesh = WingMesh(wing, 15, 30)

surf_pans = surface_panels(wing_mesh)
surf_pts  = surface_coordinates(wing_mesh)

## Freestream velocity
fs = Freestream(alpha = 5.0)
V∞ = 15.

V = V∞ * velocity(fs)

## Template:

# AIC
# [ (Nc x Ns) x (Nc x Ns) | Nw x (Nc x Ns);
#   --------------------------------------
#   Nw x (Nc x Ns)        | I_Nw          ], Note that Nw = Ns

# [ A_ff         | A_fw ;
#  ----------------------
# [ kutta ⊗ I_Nw | I_Nw ]

Nc = size(surf_pans, 1)
Ns = size(surf_pans, 2)
Nw = Ns

##

# AIC = [ rand(Ns * Nc, Ns * Nc) rand(Nc * Ns, Nw); kron(AeroFuse.DoubletSource.kutta_condition(Nc, Nw), I(Nw)) I(Nw) ]

## Aerodynamic Influence Coefficient matrix
AIC_ff = [ constant_quadrilateral_doublet_potential(1., pan_i, collocation_point(pan_j)) 
        for pan_i in surf_pans, pan_j in surf_pans ]

AIC_fw = [ constant_quadrilateral_doublet_potential(1., pan_i, p_j) 
for pan_i in surf_pans, p_j in collocation_point.(wake_panels) ]

kutta = zeros(Nw, Nc * Ns)
kutta[1:end-1,2:end] .= 1
kutta[Nc+1] ...

AIC = [ AIC_ff AIC_fw ; 
        kutta  I(Ns)  ]

# Boundary condition (no sources yet)
# BC 
# [ Nc x Ns;
#   -------
#   0_Nw   ]
RHS = -[ dot(V, p_j) for p_j in collocation_point.(surf_panels); 0 ]



# Solve
φ = AIC \ RHS

## Plotting
using Plots

plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.2, 1.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)