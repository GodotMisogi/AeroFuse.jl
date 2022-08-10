##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO

##
panel = Panel3D(
    Point3D( 1.0, -1.0,  0.0), #4
    Point3D(-1.0, -0.9,  0.0), #3
    Point3D(-1.0,  1.1,  0.0), #2
    Point3D( 1.0,  1.0,  0.0), #1
)

panel = Panel3D(
    Point3D( 0.1, -0.1,  0.0), #4
    Point3D(-0.1, -0.09,  0.0), #3
    Point3D(-0.1,  0.11,  0.0), #2
    Point3D( 0.1,  0.1,  0.0), #1
)

# point = Point3D(0.024698438432502412, -0.1209022472975474, -5.551115123125783e-17)
point = Point3D(0.2,0.2,0.2)
ϵ = 1.0e-7
pertx = Point3D(ϵ, 0., 0.)
perty = Point3D(0., ϵ, 0.)
pertz = Point3D(0., 0., ϵ)

(quadrilateral_doublet_potential(panel, point + pertx) - quadrilateral_doublet_potential(panel, point)) / ϵ
(quadrilateral_doublet_potential(panel, point + perty) - quadrilateral_doublet_potential(panel, point)) / ϵ
(quadrilateral_doublet_potential(panel, point + pertz) - quadrilateral_doublet_potential(panel, point)) / ϵ

quadrilateral_doublet_velocity(panel, point)
# quadrilateral_source_velocity_farfield(1, panel, point)

##
plot(xlabel="x", ylabel="y", zlim=(-2,3))
[plot!(p) for p in plot_panels([
    panelview[3],
    panelview[2],
    # panelview[end-2],
    # panelview[end-1],
])]
plot!(aspect_ratio=:equal)

##
plot(xlabel="x", ylabel="y")
[plot!(p) for p in plot_panels([
    surf_pans[1:3];
    surf_pans[end-2:end]
])]
plot!()

##
plot(xlabel="x", ylabel="y", zlim=(-2,3))
[plot!(p) for p in plot_panels(
    [panel, lpv3]
)]
plot!(aspect_ratio=:equal)

##
plot(xlabel="x", ylabel="y", zlim=(-2,3))
[plot!(p, color=:grey, linewidth=0.5) for p in plot_panels(
    panelview
)]
plot!(aspect_ratio=:equal, legend=:false)

##
x = LinRange(-1,1,1001)
y = (x -> quadrilateral_doublet_potential(1, lpv3, Point3D(.02,-.12,0+x))).(x)
plot(x,y)

##
av = point - coord[i];
as = norm(av);
bv = point - coord[i%4+1];
bs = norm(bv);
dv = coord[i%4+1] - coord[i];
ds = norm(dv);
al = av ⋅ l;
am = av ⋅ m;
bl = bv ⋅ l;
bm = bv ⋅ m;
dl = dv ⋅ l;
dm = dv ⋅ m;

A = dl * (am^2 + gn^2) - al * am * dm
B = dl * (bm^2 + gn^2) - bl * bm * dm

NUM = A * B + gn^2 * as * bs * dm^2
DEN = dm * gn * (bs * A - as * B)

atan(NUM/DEN)