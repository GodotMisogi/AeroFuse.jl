## Example script for linear-strength vortex panel method
using LinearAlgebra
using Base.Iterators
using Seaborn
using AeroMDAO

## Airfoil
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs     = (0., 0.)
# airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
airfoil = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = false))
V, α    = 1., 5.
ρ       = 1.225
uniform = Uniform2D(V, α)
num_pans = 80
panels  = paneller(airfoil, num_pans);

## Vortex panel
A = @time vortex_influence_matrix(panels)
b = @time neumann_boundary_condition(panels, velocity(uniform))
γs = A \ b

##
cl = 2 * sum((forward_sum(γs) / 2) .* panel_length.(panels)) / uniform.magnitude

## Velocities and pressure coefficient
# panel_vels = [ velocity(uniform) .+ sum(vortex_velocity.(γs[1:end-1], γs[2:end], panels, Ref(panel_i))) for panel_i in panels ]

qts = γs
cps = @. 1 - qts^2 / uniform.magnitude^2

## Pressure coefficient
# upper, lower = get_surface_values(panels, cps)
# upper, lower = [ upper; lower[1] ], [ lower; upper[1] ]
# plot(first.(upper), last.(upper), label = "Upper")
# plot(first.(lower), last.(lower), label = "Lower")

plot(first.(panel_points(panels)), cps)
ylim(maximum(cps), minimum(cps))
xlabel("(x/c)")
ylabel("Cp")
show()


## Plotting
x_domain, y_domain = (-1, 2), (-1, 1)
grid_size = 50
x_dom, y_dom = range(x_domain..., length = grid_size), range(y_domain..., length = grid_size)
grid = product(x_dom, y_dom)
pts = panel_points(panels);

vortex_vels = [ velocity(uniform) .+ sum(vortex_velocity.(γs[1:end-1], γs[2:end], panels, x, y)) for (x, y) in grid ]

speeds = @. sqrt(first(vortex_vels)^2 + last(vortex_vels)^2)
cps = @. 1 - speeds^2 / uniform.magnitude^2

##
contourf(first.(grid), last.(grid), cps)
CB1 = colorbar(label = "Pressure Coefficient (Cp)")
# quiver(first.(grid), last.(grid), first.(vortex_vels), last.(vortex_vels), speeds)
streamplot(first.(grid)', last.(grid)', first.(vortex_vels)', last.(vortex_vels)', color = speeds', cmap = "coolwarm", density = 2)
CB2 = colorbar(orientation = "horizontal", label = "Relative Airspeed (U/U∞)")
fill(first.(pts), last.(pts), color = "black", zorder = 3)
tight_layout()
show()
