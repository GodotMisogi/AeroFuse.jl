using AeroFuse
using Plots

## NACA
naca_foil = naca4(4,6,3,6) # NACA 4-digit

## Plot
plot(naca_foil, 
    aspect_ratio = 1, 
    camber = true, 
    thickness = true
)

## Fitting to Kulfan Class Shape Transformation parametrisation
num_dv = 8

# Coordinates fitting
up, low  = split_surface(naca_foil)
alpha_u  = coordinates_to_CST(up, num_dv)
alpha_l  = coordinates_to_CST(low, num_dv)
cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), (0., 0.))

## Plot
plot(cst_foil)
plot!(naca_foil)

## Control surface
cst_foil_flap = control_surface(cst_foil, 
                    angle = 5, # Angle of deflection
                    hinge = 0.75, # Normalized point of rotation along chord line
                )

##
plot!(cst_foil_flap)