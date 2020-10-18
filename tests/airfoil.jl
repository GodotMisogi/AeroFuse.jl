##
using Revise
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: read_foil, foil_camthick, camthick_foil, cosine_foil, kulfan_CST, naca4

##
foilpath = "airfoil_database/ys930.dat"

## Foil processing
coords = read_foil(foilpath)

## Cosine spacing
cos_foil = cosine_foil(coords)

## Camber-thickness transformations
xcamthick = foil_camthick(coords)
foiler = camthick_foil(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])

## CST Testing
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
alphas = [alpha_u alpha_l]
dzs = (1e-4, 1e-4)
cst_foil = kulfan_CST(alphas, dzs, 0.2)

## NACA 4-digit airfoils
naca = naca4((2,4,1,2))

## Plotting library
using Plots
plotly();

## Cosine and camber-thickness
plot(cos_foil[:,1], cos_foil[:,2],
     label = "Cosine")
plot!(xcamthick[:,1], xcamthick[:,2],
     label = "Camber")
plot!(xcamthick[:,1], xcamthick[:,3], 
     label = "Thickness", aspect_ratio = :equal)


## CST
plot(cst_foil[:,1], cst_foil[:,2], 
    label = "CST", aspect_ratio = :equal)

## NACA
plot(naca[:,1], naca[:,2], 
    label = "NACA", aspect_ratio = :equal)