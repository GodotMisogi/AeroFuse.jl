using AeroFuse
using Plots

## NACA
naca_foil = naca4(2,4,1,2) # NACA 4-digit


## 
plot(naca_foil)

## Fitting to Kulfan Class Shape Transformation parametrisation
num_dv = 8

# Coordinates fitting
up, low  = split_surface(naca_foil)
alpha_u  = coordinates_to_CST(up, num_dv)
alpha_l  = coordinates_to_CST(low, num_dv)
cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), (0., 0.))

##
plot(cst_foil)