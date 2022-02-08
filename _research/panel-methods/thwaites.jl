xs = range(0, 10, length = 100)
dxs = xs[2:end] .- xs[1:end-1]
ue = 1 .- exp.(-xs)
due_dx = exp.(-xs)

ν = 1E-5;
λ0 = 0.225 / 2.65;
θ2 = λ0 * ν / due_dx[1];

θs = [sqrt(θ2)];
for i = 2:length(xs)-1
    rhs = θ2 / dxs[i] + 0.225 * 2 * ν / ue[i];
    θ2 = rhs / (1 / dxs[i] + 2 * 2.65 / ue[i] * due_dx[i]);
    θs = [θs; sqrt(θ2)];
end
θs

##
using Plots
plot(xs[1:end-1], θs, xlabel = "x", ylabel = "δ*")


##
xs      = range(0, 10, length = 100)
dxs     = xs[2:end] .- xs[1:end-1]
ue      = 1 .- exp.(-xs)
due_dx  = exp.(-xs)

ν       = 1e-5;
λ0      = 0.225 / 2.65;
θ2      = λ0 * ν / due_dx[1];

δ_thick(θ2, dx, ue, ν) = θ2 / dx + 0.225 * 2 * ν / ue;
θ_thick(δ, dx, ue, due_dx) = δ / (1 / dx + 2 * 2.65 / ue * due_dx);
eval_θ(θ2, dx, ue, due_dx, ν) = θ_thick(δ_thick(θ2, dx, ue, ν), dx, ue, due_dx)
eval_θ(θ2, (dx, ue, due_dx, ν)) = eval_θ(θ2, dx, ue, due_dx, ν)


##
iterate(f, x0, ps, n) = foldl(∘, fill(x -> f(x, ps), n))(x0)
iterator(f, x0, ps, n) = iterate.(Ref(f), Ref(x0), ps[1:n], 1:n)

##
θs = iterator(eval_θ, θ2, (collect ∘ zip)(dxs[2:end], ue[2:end-1], due_dx[2:end-1], fill(ν, length(dxs[2:end]))), 98)

sqrt.(θs)

##
using Plots
plot(xs[1:end-2], sqrt.(θs), xlabel = "x", ylabel = "δ*", label = :none)
