include("../src/AeroMDAO.jl")
include("../src/FoilParametrization.jl")

using .AeroMDAO: Point2D, Point3D, Panel
using .FoilParametrization: CSTBase, shape_function, cst_coords, naca4
using PyPlot

# 2D Panel test
ps = Point2D.(0:5, 1:6)

panels = Panel(ps)

# 3D Panel test
ps = Point3D.(0:5, 1:6, 2:7)

panels = Panel(ps)

# NACA 4-digit airfoils
coords = naca4((2,4,1,2))

plot(coords[:,1], coords[:,2])
axis("equal")
show()
# CST test
