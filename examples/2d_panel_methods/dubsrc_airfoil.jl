## Example script for constant-strength doublet-source panel method
using AeroMDAO

##
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs     = (0., 0.)
airfoil = kulfan_CST(alpha_u, alpha_l, dzs, (0., 0.), 60);      # Kulfan CST airfoil
# airfoil = naca4(((0,0,1,2), 100; sharp_trailing_edge = true)) # NACA 4-digit airfoil
uniform = Uniform2D(1., 0.)
system  = @time solve_case(
                     airfoil, uniform;
                     num_panels = 80
                    );

##
panels    = system.surface_panels
@time ues = surface_velocities(system);
@time cl  = lift_coefficient(system)
@time cls, cms, cps = surface_coefficients(system)

## Printing
println("Cl: $cl")
println("Σᵢ Clᵢ: $(sum(cls))")
println("Σᵢ Cmᵢ: $(sum(cms))")

## Plotting libraries
using Plots
plotlyjs();

## Pressure coefficients
plot((first ∘ collocation_point).(panels[1:end-1]), cps, marker = 2, label = :none, yflip = true, xlabel = "(x/c)", ylabel = "Cp")

## Lift distribution
plot((first ∘ collocation_point).(panels[1:end-1]), cls, marker = 2, label = :none, xlabel = "(x/c)", ylabel = "Cl")

## Angle of attack sweep
function alpha_sweep(α, airfoil)
    uniform = Uniform2D(1.0, α)
    solve_case(airfoil, uniform, num_panels = 80)[1]
end

αs = 0:10
cls = @time alpha_sweep.(αs, Ref(airfoil))