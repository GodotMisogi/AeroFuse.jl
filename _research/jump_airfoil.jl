## 
using Revise
using AeroMDAO
using JuMP
using Ipopt
using Plots
# plotlyjs()

cl0 = 0.5      # Target lift coefficient

## Lower surface optimization
#============================================#

# Objective function
function optimize_lower(α, α_l...)
     α_u = [0.2, 0.3, 0.2, 0.15, 0.2]
     dzs = (0., 0.)
     airfoil = (Foil ∘ kulfan_CST)(α_u, [ α_l... ], dzs, 0., 80)
     uniform = Uniform2D(1.0, α)
     cl = solve_case(airfoil, uniform, num_panels = 60)[1]
end

## Test
α_u = [0.2, 0.3, 0.2, 0.15, 0.2]
dzs = (0., 0.)
α_l0 = [-0.2, -0.1, -0.1, -0.001]
lower_cl = optimize_lower(3.0, α_l0...)
println("Lower Cl: $lower_cl")

## Optimization
lower_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 5

@variable(lower_design, -1. <= α_l[1:num_dv] <= 0.) 
@variable(lower_design, -5. <= α <= 10.) 

register(lower_design, :optimize_lower, num_dv + 1, optimize_lower, autodiff = true)

@NLobjective(lower_design, Min, abs(optimize_lower(α, α_l...) - cl0))

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


## Camber optimization
#============================================#

# Parameters
α_t = [0.26, 0.23, 0.43, 0.13]
dct = (0., 0.)

# Camber optimization
function optimize_camber(α, α_c...)
    airfoil = (Foil ∘ camber_CST)([ α_c... ], α_t, dct, 0., 80)
    uniform = Uniform2D(1.0, α)
    cl = solve_case(airfoil, uniform, 60)[1]
end

## Test
α_c0 = [0.0066, 0.04, -0.06, 0.25]
camber_cl = optimize_camber(3.0, α_c0...)
println("Camber Cl: $camber_cl")

## Optimization model
cam_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 4

@variable(cam_design, -0.1 <= α_c[1:num_dv] <= 0.1) 
@variable(cam_design, -5. <= α <= 10.) 

register(cam_design, :optimize_camber, num_dv + 1, optimize_camber, autodiff = true)

@NLobjective(cam_design, Min, abs(optimize_camber(α, α_c...) - cl0))

## Run optimization
optimize!(cam_design)

## Print
println("Angle: $(value(α)), Optical Camber Cl: $(objective_value(cam_design))")
α_c_vals = value.(α_c)

## Plot
camber_test = camber_CST(α_c0, α_t, dct, 0., 80)
plot(camber_test[:,1], camber_test[:,2], 
     label = "Test Camber", aspectratio = 1)

camber_opt = camber_CST([ α_c_vals[1:end-1]... ], α_t, dct, 0., 80)
plot!(camber_opt[:,1], camber_opt[:,2], 
     label = "Optimal Camber", aspectratio = 1)


## NACA 4-series optimization
#============================================#

# Objective function
function optimize_naca4(α, m, p, t, c)
    airfoil = naca4((m, p, t, c), 40, sharp_trailing_edge = true)
    panels = make_2Dpanels(airfoil)

    uniform = Uniform2D(1.0, α)
    cl = solve_case(panels, uniform)
end

## Test
coeffs_0 = (2, 4, 1, 2)
naca_cl = optimize_naca4(3.0, coeffs_0...)
println("NACA4 Cl: $naca_cl")

## Optimization model
naca_design = Model(with_optimizer(Ipopt.Optimizer))

@variable(naca_design, 1e-4 <= digits[1:4] <= 2.) 
@variable(naca_design, -5. <= α <= 15.) 

register(naca_design, :optimize_naca4, 5, optimize_naca4, autodiff = true)

@NLobjective(naca_design, Max, abs(optimize_naca4(α, digits...) - cl0))

## Run optimization
optimize!(naca_design)

## Print
println("Angle: $(value(α)), Optimal NACA Cl: $(objective_value(naca_design))")
naca_vals = value.(digits)

## Plot
naca_test = naca4(coeffs_0, 40, sharp_trailing_edge = true)
plot(naca_test[:,1], naca_test[:,2], 
     label = "Test NACA", aspectratio = 1)

naca_opt = naca4(tuple(naca_vals...), 40, sharp_trailing_edge = true)
plot!(naca_opt[:,1], naca_opt[:,2], 
     label = "Optimal NACA", aspectratio = 1)
