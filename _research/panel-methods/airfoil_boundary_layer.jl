##
using AeroMDAO
using LinearAlgebra

##
alpha_u  = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l  = [-0.2, -0.1, -0.1, -0.001]
dzs      = (0.,0.)
# airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
airfoil  = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = true))
# airfoil = Foil(pts)
V, α     = 6., 0.
uniform  = Uniform2D(V, α)
num_pans = 64
num_wake = 28
u        = velocity(uniform)

## Doublet-source panel method evaluation
@time cl, cls, cms, cps, panels =
solve_case(airfoil,
           uniform;
           viscous     = false,
           sources     = false,
           wake_length = 1e5,
           wake_panels = num_wake,
           num_panels  = num_pans)

@show cl
@show sum(cls)
@show sum(cms)
@show re = reynolds_number(1.225, V, 1.0, 1.7894e-5)

## Panel method setup
airfoil = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = true))
panels  = paneller(airfoil, num_pans);
wake    = wake_panel(panels, 1e5, deg2rad(α))
wakes   = wake_panels(panels, 1.0, 1., num_wake)

## Evaluate inviscid flow
inv_res  = solve_inviscid_vortices(panels, wakes, velocity(uniform))

## Evaluate viscous flow
visc_res = solve_viscous_case(panels, wakes, uniform)

##

##
using Seaborn
# plot(getindex.(collocation_point.(panels), 1), cps, yflip = true)

pts = panel_points(panels)
wake_pts = panel_points(wakes)
plot(getindex.(pts, 1), getindex.(pts, 2), marker = :dot, markersize = 2, aspect_ratio = 1)
plot(getindex.(wake_pts, 1), getindex.(wake_pts, 2), marker = :dot, markersize = 2)
show()