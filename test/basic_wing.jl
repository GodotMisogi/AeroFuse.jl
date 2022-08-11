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
wing_mesh = WingMesh(wing, [20], 10; span_spacing = Uniform())

surf_pts  = surface_coordinates(wing_mesh)
left_edge = (surf_pts[:,1] + reverse(surf_pts[:,1])) / 2
right_edge = (surf_pts[:,end] + reverse(surf_pts[:,end])) / 2
surf_pts = [left_edge surf_pts right_edge]

surf_pans = make_panels(surf_pts)

# Freestream velocity
α = 8.0; β = 0.0; Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
@time prob = solve_system_neumann(surf_pans, Umag, fs, 50)
φs = permutedims(
    reshape(
        prob.singularities[1:length(prob.surface_panels)], 
        size(prob.surface_panels,2), size(prob.surface_panels,1))
)

##
vs = surface_velocities(prob)
@time cls, cps = surface_coefficients(prob, wing);
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots
plotly()
plt_surfs = plot_panels(surf_pans)

xl=-1;xh=4;yl=-3;yh=3;zl=-1;zh=1;
p = Plots.plot(; xlabel="x",ylabel="y",zlabel="z", grid=:true, xlim = (xl, xh), ylim = (yl, yh), zlim=(zl, zh))

for pan in plt_surfs
    p = Plots.plot!(pan, color = :grey)
end

plt_surfs = plot_panels(prob.wake_panels)
for pan in plt_surfs
    p = Plots.plot!(pan, color = :blue)
end

# plot!(
#     getindex.(collocation_point.(surf_pans), 1),
#     getindex.(collocation_point.(surf_pans), 2),
#     getindex.(collocation_point.(surf_pans), 3),
#     st = :surface
# )

plot!(
    getindex.(surf_pts, 1),
    getindex.(surf_pts, 2),
    getindex.(surf_pts, 3),
    st = :surface)

plot!(p)
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

##
