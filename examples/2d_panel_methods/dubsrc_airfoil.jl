## Example script for constant-strength doublet-source panel method
using AeroMDAO

##
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs     = (0., 0.)
# airfoil = kulfan_CST(alpha_u, alpha_l, dzs, (0., 0.), 60);      # Kulfan CST airfoil
airfoil = naca4((0,0,1,2), 100; sharp_trailing_edge = true) # NACA 4-digit airfoil
uniform = Uniform2D(angle = 5)
prob  = @time solve_case(
                     airfoil, uniform;
                     num_panels = 80
                    );

##
panels    = prob.surface_panels
@time ues = surface_velocities(prob);
@time cl  = lift_coefficient(prob)
@time cls, cms, cps = surface_coefficients(prob);

## Printing
println("Cl: $cl")
println("Σᵢ Clᵢ: $(sum(cls))")
println("Σᵢ Cmᵢ: $(sum(cms))")

## Plotting libraries
using Plots
gr();

## Pressure coefficients

# Split surfaces and values
cp_upper, cp_lower = get_surface_values(panels, cps)

plot(marker = 2, label = :none, yflip = true, xlabel = "(x/c)", ylabel = "Cp")
plot!(cp_upper, label = "Upper")
plot!(cp_lower, label = "Lower")

## Lift distribution
cl_upper, cl_lower = get_surface_values(panels, cls)

plot(marker = 2, label = :none, yflip = true, xlabel = "(x/c)", ylabel = "Cl")
plot!(cl_upper, label = "Upper")
plot!(cl_lower, label = "Lower")

## Angle of attack sweep example
function alpha_sweep(α, airfoil)
    uniform = Uniform2D(1.0, α)
    lift_coefficient(solve_case(airfoil, uniform, num_panels = 80))
end

αs = 0:10
cls = @time alpha_sweep.(αs, Ref(airfoil))

##
plot(αs, cls, ylabel = "Cl", xlabel = "α (deg)", label = :none)