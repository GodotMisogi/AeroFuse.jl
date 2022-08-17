## Example script for linear-strength vortex panel method
using LinearAlgebra
using Base.Iterators
using AeroFuse

## Airfoil
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs     = (0., 0.)
# airfoil = kulfan_CST(alpha_u, alpha_l, dzs, 0.0, 60);
airfoil = naca4((0,0,1,2))
V, α    = 1., 5.
ρ       = 1.225
uniform = Uniform2D(V, α)
num_pans = 80
panels  = make_panels(airfoil, num_pans);

## Vortex panel
A = @time vortex_influence_matrix(panels)
b = @time neumann_boundary_condition(panels, velocity(uniform))
γs = A \ b

##
cl = 2 * sum(@. (γs[1:end-1] + γs[2:end]) / 2 * panel_length(panels)) / uniform.magnitude

## Velocities and pressure coefficient
# panel_vels = [ velocity(uniform) .+ sum(vortex_velocity.(γs[1:end-1], γs[2:end], panels, Ref(panel_i))) for panel_i in panels ]

## Plotting
using CairoMakie

CairoMakie.activate!()

using LaTeXStrings

const LS = LaTeXString

## Pressure coefficient
qts = γs
cps = @. 1 - qts^2 / uniform.magnitude^2

upper, lower = get_surface_values(panels, cps)
upper, lower = [ lower[end]; upper ], [ upper[end]; lower ]

# Plot
f1 = Figure(resolution = (900, 400))
ax = Axis(f1[1,1], aspect = 1, yreversed = true, xlabel = L"(x/c)", ylabel = L"C_p")
l1 = lines!(upper, label = L"\mathrm{Upper}")
l2 = lines!(lower, label = L"\mathrm{Lower}")

axislegend(L"\mathrm{Surface}")

# Velocity field
x_domain, y_domain = (-1, 2), (-1, 1)
grid_size = 50
x_dom, y_dom = range(x_domain..., length = grid_size), range(y_domain..., length = grid_size)
grid = product(x_dom, y_dom)
pts = panel_points(panels);

total_velocity(xs, ys) = map((x, y) -> Point2f(velocity(uniform) .+ sum(vortex_velocity.(γs[1:end-1], γs[2:end], panels, x, y))...), xs, ys)

vortex_vels = total_velocity(first.(grid), last.(grid))
speeds = @. sqrt(first(vortex_vels)^2 + last(vortex_vels)^2)
cps = @. 1 - speeds^2 / uniform.magnitude^2

# Plot
ax2 = Axis(f1[1,2], aspect = 1)
hm = contourf!(x_dom, y_dom, speeds, levels = 50)

streamplot!(total_velocity, first.(grid), last.(grid), colormap = Reverse(:viridis))
poly!(pts, color = "black")

Colorbar(f1[1,3], hm, label = L"Normalized Speed $(U/U_\infty)$")
colsize!(f1.layout, 2, Aspect(1, 1.0))

# f1[0,:] = Label(f1, LS("Linear Vortex Panel Method"))

# arrows!(x_dom, y_dom, first.(vortex_vels), last.(vortex_vels))
# CB2 = colorbar(orientation = "horizontal", label = "Relative Airspeed (U/U∞)")
# tight_layout()
# show()
f1

##
save("plots/LinearVortex.svg", f1, px_per_unit = 1.5)