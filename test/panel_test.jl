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