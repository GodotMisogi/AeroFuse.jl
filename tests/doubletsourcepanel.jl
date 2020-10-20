##
using Revise
includet("../src/PanelMethods.jl")
includet("../src/FoilParametrization.jl")
includet("../src/MathTools.jl")

##
using .PanelMethods: make_panels, collocation_point, split_panels, panels_xs, panels_ys, solve_case, Uniform2D, grid_data, pressure_coefficient
using .FoilParametrization: read_foil, cosine_foil, kulfan_CST, naca4
using .MathTools: linspace, ×, <<
using BenchmarkTools

##
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001, -0.02]
alphas = [alpha_u alpha_l]
dzs = (1e-4, 1e-4)
airfoil = kulfan_CST(alphas, dzs, 0.0, 100)

##
uniform = Uniform2D(1.0, 5.0)
panels = make_panels(airfoil)
lower_panels, upper_panels = split_panels(panels)

##
@time dub_src_panels, cl = solve_case(panels, uniform)
# println("Lift Coefficient: $cl")

##
cls = []
for alpha in 0:15
    uniform = Uniform2D(5.0, alpha)

    @time dub_src_panels, cl = solve_case(panels, uniform)
    println("Angle: $alpha, Lift Coefficient: $cl")

    append!(cls, cl)
end

## Plotting libraries
using Plots, LaTeXStrings
plotlyjs();

## Lift polar
plot(0:15, cls, 
        xlabel = L"\alpha", ylabel = L"C_l")

## Plotting domain
x_domain, y_domain = (-1, 2), (-1, 1)
grid_size = 50
x_dom, y_dom = linspace(x_domain..., grid_size), linspace(y_domain..., grid_size)
grid = x_dom × y_dom

vels, pots = grid_data(dub_src_panels, grid)
cp = pressure_coefficient.(uniform.mag, vels);

lower_panels, upper_panels = split_panels(dub_src_panels);

## Airfoil plot
plot( panels_xs(upper_panels), panels_ys(upper_panels), 
        label = "Upper", markershape = :circle,
        xlabel = "x", ylabel = "C_p")
plot!(panels_xs(lower_panels), panels_ys(lower_panels),
        label = "Lower", markershape = :circle,
        xlabel = "x", ylabel = "C_p")

## Pressure coefficient
plot( panels_xs(upper_panels), :cp .<< upper_panels, 
        label = "Upper", markershape = :circle, 
        xlabel = "x", ylabel = "C_p")
plot!(panels_xs(upper_panels), :cp .<< lower_panels, 
        label = "Lower", markershape = :circle, yaxis = :flip)

## Control volume
p1 = contour(x_dom, y_dom, cp, fill = true)
plot(p1)
plot!(first.(:start .<< panels), last.(:start .<< panels), 
      color = "black", label = "Airfoil", aspect_ratio = :equal)