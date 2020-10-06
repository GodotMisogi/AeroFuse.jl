include("../src/FoilParametrization.jl")

using .FoilParametrization: Foil, read_foil, linspace, foil_camthick, camthick_foil, cosine_foil, kulfan_CST
using PyPlot

foilpath = "airfoil_database/ys930.dat"

# Foil processing
coords = read_foil(foilpath)

# Cosine spacing
cos_foil = cosine_foil(coords)

# Camber-thickness transformations
xcamthick = foil_camthick(coords)
foiler = camthick_foil(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])

# CST Testing
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
alphas = [alpha_u alpha_l]
dzs = (1e-4, 1e-4)
cst_foil = kulfan_CST(alphas, dzs, 0.2)

# Plotting
figure(1)
plot(cos_foil[:,1], cos_foil[:,2])
plot(xcamthick[:,1], xcamthick[:,2], marker=".")
plot(xcamthick[:,1], xcamthick[:,3], marker="x")
axis("equal")

figure(2)
plot(cst_foil[:,1], cst_foil[:,2])
axis("equal")
show()
