##
using Revise
using AeroMDAO
using JuMP
using Ipopt
using Plots
using Xfoil: set_coordinates, solve_alpha
plotlyjs()

cl0 = 0.5      # Target lift coefficient

## Lower surface optimization
#============================================#

# Objective function
function coords(α_l...)
    α_u = [0.1, 0.1, 0.1, 0.1, 0.1]
    dzs = (0., 0.)
    airfoil = (Foil ∘ kulfan_CST)(α_u, [ α_l... ], dzs, 0., 80).coords
end

function solve(α, α_l...)
    airfoil = coords(α_l...)
    set_coordinates(airfoil[:,1], airfoil[:,2])
    solve_alpha(α, reynolds_number(1.225, 1., 1., 1.7894e-5))
end

cd(α, α_l...) = solve(α, α_l...)[1]
cl(α, α_l...) = solve(α, α_l...)[2]

## Test
dzs = (0., 0.)
α_u0 = [0.2, 0.3, 0.2, 0.15, 0.2]
α_l0 = [-0.1, -0.1, -0.1, -0.001]
lower = solve(0.0, α_l0...)
println("Lower: $lower")

## Optimization
lower_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 5

@variable(lower_design, -1. <= α_l[1:num_dv] <= 0.)
@variable(lower_design, -5. <= α <= 10.)

register(lower_design, :optimize_cd, num_dv + 1, cd, autodiff = true)
register(lower_design, :cons_lift, num_dv + 1, cd, autodiff = true)

@NLobjective(lower_design, Min, optimize_cd(α, α_l...))
@NLconstraint(lower_design, cons_lift(α, α_l...) == 0.5)

## Run optimization
optimize!(lower_design)

## Print
println("Angle: $(value(α)), Optimal Lower Cl: $(objective_value(lower_design))")
α_l_vals = value.(α_l)

## Plot
lower_test = kulfan_CST(α_u, α_l0, dzs, 0., 80)
plot(lower_test[:,1], lower_test[:,2],
     label = "Test Lower Surface", aspectratio = 1)

lower_opt = kulfan_CST(α_u, [ α_l_vals[1:end-1]... ], dzs, 0., 80)
plot!(lower_opt[:,1], lower_opt[:,2],
     label = "Optimal Lower Surface", aspectratio = 1)