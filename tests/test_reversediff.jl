## 
using Revise
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: naca4, kulfan_CST
using AeroMDAO

## 2D doublet-source panel method
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)

airfoil = kulfan_CST(alpha_u, alpha_l, (0., 0.), 0., 80)
panels = make_2Dpanels(airfoil)

##
function optimize_CST(alpha_u, alpha_l)
    
    airfoil = kulfan_CST(alpha_u, alpha_l, (0., 0.), 0., 80)
    panels = make_2Dpanels(airfoil)

    uniform = Uniform2D(1.0, 5.0)
    cl = solve_case(panels, uniform)
end

optimize_CST(alpha_u, alpha_l)

##
using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults

##
# Crashes due to Inexact Error Int(0.5)
const f_tape = GradientTape(optimize_CST, (alpha_u, alpha_l))
const compiled_f_tape = compile(f_tape)

##
inputs = (alpha_u, alpha_l)
results = (similar(alpha_u), similar(alpha_l))
all_results = map(DiffResults.GradientResult, results)
cfg = GradientConfig(inputs)

##
gradient!(results, compiled_f_tape, inputs)