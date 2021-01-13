##
using Revise
using BenchmarkTools
using TimerOutputs
using AeroMDAO

##
reset_timer!();

alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs     = (1e-4, 1e-4)
airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
uniform = Uniform2D(1., 5.)
@time cl = solve_case(airfoil, uniform, 60)

#
println("Lift Coefficient: $cl")

print_timer();

##
function alpha_sweep(α, airfoil)
    uniform = Uniform2D(1.0, α)
    solve_case(airfoil, uniform, 60)
end

##
αs = 0:10
@benchmark cls = alpha_sweep.(αs, Ref(airfoil))

## Overall

function optimize_CST(alpha_u, alpha_l)
    @timeit "Make Airfoil" airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0., 80)
    @timeit "Make Uniform2D" uniform = Uniform2D(1.0, 5.0)
    @timeit "Solve Case" solve_case(airfoil, uniform, 60)
end

##
reset_timer!();

@time optimize_CST(alpha_u, alpha_l)

print_timer();

## Plotting libraries
using Plots
plotlyjs();

## Lift polar
plot(αs, cls, xlabel = "α", ylabel = "CL")