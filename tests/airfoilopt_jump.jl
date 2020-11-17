## 
using Revise

##
includet("../src/FoilParametrization.jl")
using .FoilParametrization: naca4, kulfan_CST, camber_CST
using AeroMDAO

# Optimization libraries
using JuMP
using Ipopt

# Plotting
using Plots
plotlyjs()

##----------------------Lower surface optimization----------------------##

# Parameters
α_u = [0.2, 0.3, 0.2, 0.15, 0.2]
dzs = (0., 0.)

# Objective function
function optimize_lower(α, α_l...)
    airfoil = kulfan_CST(α_u, [ α_l... ], dzs, 0., 80)
    panels = make_2Dpanels(airfoil)

    uniform = Uniform2D(1.0, α)
    cl = solve_case(panels, uniform)
end

## Test
α_l = [-0.2, -0.1, -0.1, -0.001]
lower_cl = optimize_lower(3.0, α_l...)
println("Lower Cl: $lower_cl")

## Optimization
lower_design = Model(Ipopt.Optimizer)
num_dv = 5

@variable(lower_design, -1. <= α_l[1:num_dv] <= 0.) 
@variable(lower_design, -5. <= α <= 15.) 

register(lower_design, :optimize_lower, num_dv + 1, optimize_lower, autodiff = true)

@NLobjective(lower_design, Max, optimize_lower(α, α_l...))

## Run optimization
optimize!(lower_design)

## Print
println("Optimal Lower Cl: $objective_value(lower_design)")
α_l_vals = value.(α_l)

## Plot
lower_opt = kulfan_CST(α_u, [ α_l_vals[1:end-1]... ], dzs, 0., 80)
plot(lower_opt[:,1], lower_opt[:,2], 
     label = "Optimal Lower Surface", aspectratio = 1)


##--------------------Camber distribution optimization-------------------------#

# Parameters
α_t = [0.26, 0.23, 0.43, 0.13]
dct = (0., 0.)

# Camber optimization
function optimize_camber(α, α_c...)
    airfoil = camber_CST([ α_c... ], α_t, dct, 0., 80)
    panels = make_2Dpanels(airfoil)

    uniform = Uniform2D(1.0, α)
    cl = solve_case(panels, uniform)
end

## Test
α_c = [0.0066, 0.04, -0.06, 0.25]
camber_cl = optimize_camber(3.0, α_c...)
println("Camber Cl: $camber_cl")

## Optimization model
cam_design = Model(Ipopt.Optimizer)
num_dv = 8

@variable(cam_design, 0. <= α_c[1:num_dv] <= 0.3) 
@variable(cam_design, -5. <= α <= 15.) 

register(cam_design, :optimize_camber, num_dv + 1, optimize_camber, autodiff = true)

@NLobjective(cam_design, Max, optimize_camber(α, α_c...))

## Run optimization
optimize!(cam_design)

## Print
println("Optimcal Camber Cl: $objective_value(cam_design)")
α_c_vals = value.(α_c)

## Plot
camber_opt = camber_CST([ α_c_vals[1:end-1]... ], α_t, dct, 0., 80)
plot(camber_opt[:,1], camber_opt[:,2], 
     label = "Camber", aspectratio = 1)


##------------------------NACA 4-ditit optimization----------------------------------##

# Objective function
function optimize_naca4(α, m, p, t, c)
    airfoil = naca4((m, p, t, c), 40, sharp_trailing_edge = true)
    panels = make_2Dpanels(airfoil)

    uniform = Uniform2D(1.0, α)
    cl = solve_case(panels, uniform)
end

## Test
naca_test = (2, 4, 1, 2)
naca_cl = optimize_naca4(3.0, naca_test...)
println("NACA4 Cl: $naca_cl")

## Optimization model
naca_design = Model(Ipopt.Optimizer)

@variable(naca_design, 1e-4 <= digits[1:2] <= 2.) 
@variable(naca_design, -5. <= α <= 15.) 

register(naca_design, :optimize_naca4, 5, optimize_naca4, autodiff = true)

@NLobjective(naca_design, Max, optimize_naca4(α, 2, 4, digits...))

## Run optimization
optimize!(naca_design)

## Print
println("Optimal NACA Cl: $objective_value(naca_design)")
naca_vals = [ value.(digits); value(α) ]

## Plot
naca = naca4((2, 4, naca_vals[1], naca_vals[2]), 40, sharp_trailing_edge = true)
plot(naca[:,1], naca[:,2], 
     label = "NACA", aspectratio = 1)
