using DifferentialEquations
using Plots

const g = 9.81
L = 1.0
tspan = (0.0,pi/2)

function simplependulum!(du,u,p,t)
    θ  = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g/L)*sin(θ)
end

function bc1!(residual, u, p, t)
    residual[1] = u[end÷2][1] + pi/2 # the solution at the middle of the time span should be -pi/2
    residual[2] = u[end][1] - pi/2 # the solution at the end of the time span should be pi/2
end

bvp1 = BVProblem(simplependulum!, bc1!, [pi/2,pi/2], tspan)
sol1 = solve(bvp1, GeneralMIRK4(), dt=0.05)

plot(sol1)