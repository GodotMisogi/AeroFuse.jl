##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO

##
panel = Panel3D(
    Point3D( 1.0, -1.0,  0.0), #4
    Point3D(-1.0, -1.0,  0.0), #3
    Point3D(-1.0,  1.0,  0.0), #2
    Point3D( 1.0,  1.0,  0.0), #1
)

panel1 = Panel3D(
    Point3D( 1.0, -1.0,  1.0), #4
    Point3D(-1.0, -1.0,  1.0), #3
    Point3D(-1.0,  1.0,  1.0), #2
    Point3D( 1.0,  1.0,  1.0), #1
)

point1 = Point3D(0., 0., 10.)

ϵ = 1.0e-5
pertx = Point3D(ϵ, 0., 0.)
perty = Point3D(0., ϵ, 0.)
pertz = Point3D(0., 0., ϵ)

(quadrilateral_doublet_potential(1, panel, point + pertx) - quadrilateral_doublet_potential(1, panel, point)) / ϵ
(quadrilateral_doublet_potential(1, panel, point + perty) - quadrilateral_doublet_potential(1, panel, point)) / ϵ
(quadrilateral_source_potential(1, panel, point + pertz) - quadrilateral_source_potential(1, panel, point)) / ϵ

quadrilateral_source_velocity(1, panel, point)
quadrilateral_source_velocity_farfield(1, panel, point)

quadrilateral_source_velocity(1, panel, Point3D(0,0,1e-7))

quadrilateral_doublet_potential(1, panel, Point3D(2,0,0))

##

plot!(xs(foil_pans[npansp+1]), ys(foil_pans[npansp+1]), zs(foil_pans[npansp+1]), dpi=300, xlabel="x", ylabel="y", zlabel="z")
plot!(xs(foil_pans[npansp+2]), ys(foil_pans[npansp+2]), zs(foil_pans[npansp+2]), dpi=300, xlabel="x", ylabel="y", zlabel="z")
plot!(xs(foil_pans[npansp+3]), ys(foil_pans[npansp+3]), zs(foil_pans[npansp+3]), dpi=300, xlabel="x", ylabel="y", zlabel="z")

##
plot(xs(thispan), ys(thispan), zs(thispan))
scatter!([thispoint.x],[thispoint.y],[thispoint.z])

