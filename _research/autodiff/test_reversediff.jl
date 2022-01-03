##
using ReverseDiff: JacobianTape, JacobianConfig, GradientTape, GradientConfig, jacobian, jacobian!, gradient, gradient!, compile, DiffResults
using AeroMDAO

## 2D doublet-source panel method

alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)

## Objective function
function optimize_CST(alpha_u, alpha_l)
    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0., 80)
    uniform = Uniform2D(1.0, 5.0)
    cl = solve_case(airfoil, uniform, 60)
end

optimize_CST(alpha_u, alpha_l)

## Gradient
const CST_grad_tape = GradientTape(optimize_CST, (alpha_u, alpha_l))
const compiled_CST_grad_tape = compile(CST_grad_tape)

##
inputs = (alpha_u, alpha_l)
grad_results = similar.(inputs)
grad_all_results = map(DiffResults.GradientResult, grad_results)
cfg = GradientConfig(inputs)

gradient!(grad_results, compiled_CST_grad_tape, inputs)

## Jacobian
const CST_jac_tape = JacobianTape(optimize_CST, (alpha_u, alpha_l))
const compiled_CST_jac_tape = compile(CST_jac_tape)

##
inputs = ([0.2, 0.3, 0.2, 0.15, 0.2], [-0.2, -0.1, -0.1, -0.001])
jac_results = similar.(inputs)
jac_all_results = map(DiffResults.JacobianResult, jac_results)
cfg = JacobianConfig(inputs)

jacobian!(jac_results, compiled_CST_jac_tape, inputs)