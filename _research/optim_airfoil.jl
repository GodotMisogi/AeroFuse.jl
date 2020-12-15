## 
using Revise
using AeroMDAO
using Optim
using Plots
plotlyjs()

## Surface optimization
#============================================#

# Objective function
function optimize_CST(αs, α, n_upper :: Integer)
    dzs = (0., 0.)
    airfoil = (Foil ∘ kulfan_CST)(αs[1:n_upper], αs[n_upper+1:end], dzs, 0., 80)
    uniform = Uniform2D(1.0, α)
    cl = solve_case(airfoil, uniform, 60)
end

## Test
dzs = (0., 0.)
α_u0 = [0.1, 0.1, 0.1, 0.1, 0.1]
α_l0 = [-0.1, -0.1, -0.1, -0.1]
α0 = [α_u0; α_l0]

n_upper = length(α_u0)  # Number of upper surface variables
α = 0.                  # Angle of attack
cl0 = 0.8               # Target lift coefficient
CST_cl = optimize_CST(α0, α, length(α_u0))
println("CST Cl: $CST_cl")

## Plot
upper_test = kulfan_CST(α0[1:n_upper], α0[n_upper+1:end], dzs, 0., 80)
plot(upper_test[:,1], upper_test[:,2], 
     label = "Initial Surface", aspectratio = 1)

## Optimization
resi = optimize(x -> abs(optimize_CST(x, α, n_upper) - cl0), 
                α0, 
                method = GradientDescent(),
                autodiff = :forward,
                store_trace = true,
                show_trace = true)

## Plot
upper_opt = kulfan_CST(resi.minimizer[1:n_upper], resi.minimizer[n_upper+1:end], dzs, 0., 80)

fig = plot(aspectratio = 1, dpi = 300)
plot!(upper_test[:,1], upper_test[:,2], 
        label = "Initial Surface")
plot!(upper_opt[:,1], upper_opt[:,2], 
        label = "Optimal Surface")
# savefig(fig, "_research/tmp/cuck.pdf")
                
## Optimal test
CST_cl = optimize_CST(resi.minimizer, α, n_upper)
println("CST Cl: $CST_cl")