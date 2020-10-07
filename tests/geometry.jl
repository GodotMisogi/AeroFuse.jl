include("../src/AeroMDAO.jl")
include("../src/FoilParametrization.jl")

using .AeroMDAO: Point2D, Point3D, Panel
using .FoilParametrization: CSTBase, shape_function, cst_coords, naca4, bernstein_class, bernstein_basis, split_foil
using PyPlot

# 2D Panel test
ps = Point2D.(0:5, 1:6)

panels = Panel(ps)

# 3D Panel test
ps = Point3D.(0:5, 1:6, 2:7)

panels = Panel(ps)

# CST test
# kulfan = CSTBase(bernstein_class, bernstein_basis)
