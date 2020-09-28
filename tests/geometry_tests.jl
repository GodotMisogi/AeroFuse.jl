include("../src/AeroMDAO.jl")

using .AeroMDAO: Point2D, Point3D, Panel, midpoint

# 2D Test
ps = Point2D{Float64}.(range(0, stop=5), range(1, stop=6))

panels = Panel(ps)
println(midpoint(panels))

# 3D Test
ps = Point3D{Float64}.(range(0, stop=5), range(1, stop=6), range(2, stop = 7))

panels = Panel(ps)
println(midpoint(panels))