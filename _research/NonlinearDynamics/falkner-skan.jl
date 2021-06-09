using DifferentialEquations

##
"""
Falkner-Skan ODE:

    F = ψ/m
    U = u/uₑ
    S = τξ/(ρuₑ²δ) = τ/(mu/ξ)
    
    where m = ρₑuₑδ
"""
function falkner_skan_ODE!(du, x, p, η)
    F, U, S = x
    a       = p[1]
    du[1]   = U                                     # F' = U
    du[2]   = S                                     # U' = F'' = S
    du[3]   = -(1 + a) / 2 * F * S - a * (1 - U^2)  # S' = F''' = -F * F'' - β(1 - (F')²)

    du
end

function falkner_skan_BC!(R, x, p, η)
    R[1] = x[1][1]              # η = 0, F = 0
    R[2] = x[1][2]              # η = 0, F' = 0
    R[3] = x[end][2] - 1        # η = ηₑ, F' = 1

    R
end

##
ηs    = (0.0, 10.0)
x0    = [5.0, 2.0, 0.0] 
as    = range(-0.2, 0.2, length = 20)

##
prob  = BVProblem.(Ref(falkner_skan_ODE!), Ref(falkner_skan_BC!), Ref(x0), Ref(ηs), as)
sol   = solve.(prob, Ref(GeneralMIRK4()), dt = 2e-1)

##
using Plots
plotly(dpi = 300, size = (800, 600), legend = :left)

# using S

##
cock = Plots.plot(xaxis = "", 
                  yaxis = "η = n/δ",
                  linewidth = 2)
[ Plots.plot!(solution, label = "U/Uₑ = $a", vars = [(2, 0)]) for (solution, a) in zip(sol, as) ]
Plots.plot!()
Plots.savefig(cock, "lol.json")