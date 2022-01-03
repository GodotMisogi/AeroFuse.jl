##
using AeroMDAO

##
airfoil = Foil(naca4((0,0,1,2), 81; sharp_trailing_edge = false))
# airfoil = Foil(pts)q
V, α    = 6., 0.
ρ       = 1.225
uniform = Uniform2D(V, α)
num_pans = 12
num_wake = 28

panels  = paneller(airfoil, num_pans);
wakes   = wake_panels(panels, 1.0, 1., num_wake);

##
A = influence_matrix(panels)
b = boundary(panel_points(panels), uniform)

xs = A \ b

##
γs = xs[1:end-1]
ψ0 = xs[end];

##
lift(ρ, u, γs, cs) = ρ * u * (γs[1:end-1] .+ γs[2:end]) ./ 2 .* cs

ΔLs = lift(ρ, V, γs, panel_length.(panels))

L = sum(ΔLs)
# Cl = L / (0.5 * 1.225 * V^2)
Cl = 2 * sum(fwdsum(γs) / 2 .* panel_length.(panels)) / uniform.magnitude
println("Lift Coefficient: $Cl")

##
using Seaborn

pts = panel_points(reverse_panel.(panels)[end:-1:1])
plot(first.(pts), last.(pts))