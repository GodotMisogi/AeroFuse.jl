##
using Revise
using BenchmarkTools
using TimerOutputs
using AeroMDAO

##
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)
airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 100);

##
reset_timer!();

uniform = Uniform2D(1.0, 5.0)
@time cl = solve_case(airfoil, uniform, 60)
println("Lift Coefficient: $cl")

print_timer();

##
cls = []
for alpha in 0:15
    uniform = Uniform2D(5.0, alpha)

    @time cl = solve_case(airfoil, uniform, 60)
    println("Angle: $alpha, Lift Coefficient: $cl")

    append!(cls, cl)
end

## Plotting libraries
using Plots
plotlyjs();

## Lift polar
plot(0:15, cls, 
        xlabel = "α", ylabel = "CL")