include("../src/AeroMDAO.jl")

using .AeroMDAO

# 2D Test
ps = [ AeroMDAO.Point2D{Float64}(x, y) for (x, y) in zip(range(0, stop=5), range(1, stop=6)) ]

panels = AeroMDAO.Panel(ps)
println(AeroMDAO.midpoint(panels))

# 3D Test
ps = [ AeroMDAO.Point3D{Float64}(x, y, z) for (x, y, z) in zip(range(0, stop=5), range(1, stop=6), range(2, stop = 7)) ]

panels = AeroMDAO.Panel(ps)
println(AeroMDAO.midpoint(panels))