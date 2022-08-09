##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   tip_foil   = naca4(0,0,1,2),
                   root_chord = 1.0,
                   taper      = 1.0,
                   span       = 5.0,
                   dihedral   = 0.0,
                   sweep      = 0.0)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 5.0
β = 0.0
Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
prob = solve_system(surf_pans, Umag, fs, 100)

##
i = 1
j = 1
ps = prob.surface_panels
npancd, npansp = size(ps)
npanf = npancd * npansp
φs = permutedims(reshape(prob.singularities[1:npanf], npansp, npancd))
# l = normalize(p1(ps[i,j]) - p2(ps[i,j]))
# n = panel_normal(ps[i,j])
tr = get_transformation(ps[i,j])
nb1 = tr(collocation_point(ps[i+1,j]))
nb2 = tr(collocation_point(ps[i,j+1]))

A = Point3D(0., 0., φs[i,j])
B = Point3D(nb1[1], nb1[2], φs[i+1,j])
C = Point3D(nb2[1], nb2[2], φs[i,j+1])

N = (B - A) × (C - A)

N[3]/N[1], N[3]/N[2]


make_tuple(a, b) = (a, b)

##
vxs, vys = surface_velocities(system)
@time cls, cps = surface_coefficients(system);
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots
plt_surfs = plot_panels(surf_pans)

p = Plots.plot(; xlabel="x",ylabel="y",zlabel="z", aspect_ratio=:equal, grid=:true, zlim=(-2.0, 3.0))

for pan in wakes
    p = Plots.plot!(pan, color = :grey, aspect_ratio=:equal)
end

plot!(p)
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

##
