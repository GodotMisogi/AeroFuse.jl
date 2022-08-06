using AeroMDAO

panel = Panel3D(
	Point3D( 1.0,  1.0, 0.0),
	Point3D(-1.0,  1.0, 0.0),
	Point3D(-1.0, -1.0, 0.0),
	Point3D( 1.0, -1.0, 0.0)
)

point = Point3D(0.0, 0.0, 10.0)

thispan = surf_pans[1,1]
plot(xs(thispan), ys(thispan), zs(thispan), xs(thispan))

plot!(xs(foil_pans[npansp+1]), ys(foil_pans[npansp+1]), zs(foil_pans[npansp+1]), dpi=300, xlabel="x", ylabel="y", zlabel="z")
plot!(xs(foil_pans[npansp+2]), ys(foil_pans[npansp+2]), zs(foil_pans[npansp+2]), dpi=300, xlabel="x", ylabel="y", zlabel="z")
plot!(xs(foil_pans[npansp+3]), ys(foil_pans[npansp+3]), zs(foil_pans[npansp+3]), dpi=300, xlabel="x", ylabel="y", zlabel="z")