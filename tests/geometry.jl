##
# using Revise
# includet("../src/AeroMDAO.jl")
include("../src/Geometry.jl")

##
using .Geometry: Point2D, Point3D


##
ps = Point2D.(0:5, 1:6)
# panels = Panel(ps)

##
ps = Point3D.(0:5, 1:6, 2:7)
# panels = Panel(ps)

