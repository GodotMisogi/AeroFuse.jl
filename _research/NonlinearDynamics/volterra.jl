##
f1(x, y, α, η) = α * x - η * x * y
f2(x, y, β, λ) = -β * x + λ * x * y


function volterra!(R, xs, ps, t)
    R[1] = f1(xs[1], xs[2], ps[1,1], ps[1,2])
    R[2] = f2(xs[1], xs[2], ps[2,1], ps[2,2])

    R
end

##
α, η, β, λ = 0.1, 0.1, 0.2, 0.3
ps = [ α η ; 
       β λ ]
x0 = [ 2.0, 2.0 ]

##
using DifferentialEquations

tspan = (0.0, 500.0)
prob = ODEProblem(volterra!, x0, tspan, ps)
sol = solve(prob, maxiters = 500)

##
using Plots
plotly()

plot(sol,linewidth=2,xaxis="t",label=["x" "y"],layout=(2,1))