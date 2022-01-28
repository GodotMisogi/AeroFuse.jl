## Example script for linear-strength source panel method
using LinearAlgebra
using Base.Iterators
using Seaborn
using AeroMDAO

## Airfoil
airfoil = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = false))
V, α    = 1., 0.
ρ       = 1.225
uniform = Uniform2D(V, α)
num_pans = 80
panels  = make_panels(airfoil, num_pans);

## Linear-strength source panel
As = @time source_influence_matrix(panels)
bs = @time neumann_boundary_condition(panels, velocity(uniform))
σs = As \ bs

##
sum(σs)

## Pressure coefficient
pts = panel_points(panels)[1:end-1];
panel_vels = [ velocity(uniform) .+ sum(source_velocity.(σs[1:end-1], σs[2:end], panels, x, y)) for (x, y) in pts ]

qts = @. dot(panel_vels, panel_tangent(panels))
cps = @. 1 - qts^2 / uniform.magnitude^2

##
upper, lower = get_surface_values(panels, cps)
lower = [ upper[end]; lower; upper[1] ]

plot(first.(upper), last.(upper), label = "Upper")
plot(first.(lower), last.(lower), label = "Lower")
ylim(maximum(cps), minimum(cps))
xlabel("(x/c)")
ylabel("Cp")
legend()
show()

## Plotting
x_domain, y_domain = (-1, 2), (-1, 1)
grid_size = 50
x_dom, y_dom = range(x_domain..., length = grid_size), range(y_domain..., length = grid_size)
grid = product(x_dom, y_dom)
pts = panel_points(panels);

source_vels = [ velocity(uniform) .+ sum(source_velocity.(σs[1:end-1], σs[2:end], panels, x, y)) for (x, y) in grid ]

speeds = @. sqrt(first(source_vels)^2 + last(source_vels)^2)
cps = @. 1 - speeds^2 / uniform.magnitude^2

contourf(first.(grid), last.(grid), cps)
CB1 = colorbar(label = "Pressure Coefficient (Cp)")
# quiver(first.(grid), last.(grid), first.(source_vels), last.(source_vels), speeds)
streamplot(first.(grid)', last.(grid)', first.(source_vels)', last.(source_vels)', color = speeds', cmap = "coolwarm", density = 2.5)
CB2 = colorbar(orientation="horizontal", label = "Relative Airspeed (U/U∞)")
tight_layout()
fill(first.(pts), last.(pts), color = "black", zorder = 3)
show()
