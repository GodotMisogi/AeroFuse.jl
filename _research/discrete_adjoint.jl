using LinearAlgebra
using Zygote
using StaticArrays
using Base.Iterators

# Newton iteration
prod_size(tensor) = (prod ∘ size)(tensor)

function newton_step(u, R, ∂R∂u)
    i, j = size(R)
    l, k = size(u)
    jac = reshape(∂R∂u, (i*k, j*l))
    reshape(jac \ -R[:], size(u)) 
end

## Finite differencing
central_2nd(Tm, T, Tp) = Tm - 2T + Tp
central_2nd_diff(Tm, T, Tp, h) = central_2nd(Tm, T, Tp) / h^2
central_1st_diff(Tm, Tp, h) = (Tm - Tp) / 2h

compute_q(x, y) = 5e5 * (1 - (x - 0.5)^2) * (1 - (y - 0.5)^2)

function get_T(T_grid, i, j, T_boundary)
    if i < 1 || j < 1 || i > length(T_grid[:,1])|| j > length(T_grid[1,:])
        T_boundary
    else
        T_grid[i,j]
    end
end

stencil_1(T, i, j, T_boundary) = 
    [ get_T(T, i, j, T_boundary)   ;
      get_T(T, i-1, j, T_boundary) ;
      get_T(T, i+1, j, T_boundary) ;
      get_T(T, i, j-1, T_boundary) ;
      get_T(T, i, j+1, T_boundary) ]

# Difference operators specific to stencil_1
Δ²x(T, h) = central_2nd_diff(T[2], T[1], T[3], h)
Δ²y(T, k) = central_2nd_diff(T[4], T[1], T[5], k)
Δx(T, h) = central_1st_diff(T[3], T[2], h)
Δy(T, k) = central_1st_diff(T[5], T[4], k)

compute_k(β₁, β₂, β₃, β₄,  T) = β₁/T + β₂ + β₃ * T + β₄ * T^2
# + β₃ * T + β₄ * T^2
∇compute_k(β₁, β₂, β₃, β₄, T) = gradient(compute_k, β₁, β₂, β₃, β₄, T)

function residual_cell(T, β, q, dx, dy)
    k = compute_k(β..., first(T)) 
    ∂k∂T = (last ∘ ∇compute_k)(β..., first(T))
    
    k * (Δ²x(T, dx) + Δ²y(T, dy)) + ∂k∂T * (Δx(T, dx)^2 + Δy(T, dy)^2 ) + q
end

residual(T, β, q, dx, dy, T_boundary) = 
        [ let Ts = stencil_1(T, T_ind.I..., T_boundary);
            residual_cell(Ts, β, q[T_ind], dx, dy) end
            for T_ind in CartesianIndices(T) ]

# Derivative stuff
function solve_direct(x, u, ∂R∂x, ∂R∂u)
    ∂R∂u_sq = reshape(∂R∂u, (length(u[:]), length(u[:])))
    reshape(hcat((∂R∂u_sq \ -(∂R∂x)[:,:,i][:] for i in eachindex(x))...), (size(u)..., length(x)))
end

function solve_adjoint(u, ∂R∂u, dfdu) 
    reshape(∂R∂u, (length(u[:]), length(u[:])))' \ -(dfdu)'[:]
end

total_derivative_direct(∂f∂x, ψ, ∂f∂u) = ∂f∂x + [ sum(∂f∂u * ψ[n]) for n in eachindex(∂f∂x) ]

total_derivative_adjoint(∂f∂x, φ, ∂R∂x) = ∂f∂x + [ sum(permutedims(φ) * reshape(∂R∂x, (length(R[:]), length(∂f∂x)))[:,n]) for n in eachindex(∂f∂x) ]

∇residual_cell(Ts, βs, q, dx, dy) = gradient(residual_cell, Ts, βs, q, dx, dy)

function ∇residual_cell_T_update!(∂R∂T, grad, k, l)
    a, b, c, d = size(∂R∂T)  
    ∂R∂T[k,l,k,l] = grad[1]
    if k > 1 ∂R∂T[k,l,k-1,l] = grad[2] end
    if k < a ∂R∂T[k,l,k+1,l] = grad[3] end
    if l > 1 ∂R∂T[k,l,k,l-1] = grad[4] end
    if l < b ∂R∂T[k,l,k,l+1] = grad[5] end
end

∇residual(T, β, q, dx, dy, T_boundary) = gradient(residual, T, β, q, dx, dy, T_boundary)

∇residual_β(T, β, q₀, dx, dy, T_boundary) = 
        [ let Ts = stencil_1(T, R_ind.I..., T_boundary);
            ∇residual_cell(Ts, β, q₀[R_ind], dx, dy)[2][i] end for
            R_ind in CartesianIndices(R), i in eachindex(β) ]

function compute_residual_and_grad(T, β, q, dx, dy)
    # Residual
    R = [ let Ts = stencil_1(T, T_ind.I..., T_boundary);
            residual_cell(Ts, β, q[T_ind], dx, dy) end
            for T_ind in CartesianIndices(T) ]

    # Jacobian
    T_size = size(T)
    ∂R∂T = zeros(T_size..., T_size...)
    for T_ind in CartesianIndices(T)
        Ts = stencil_1(T, T_ind.I..., T_boundary)
        grad = ∇residual_cell(Ts, β, q₀[T_ind], dx, dy)[1]
        ∇residual_cell_T_update!(∂R∂T, grad, T_ind.I...)
    end

    R, ∂R∂T
end


## Grid setup
dx = 0.1
dy = 0.2
xs, ys = -1:dx:1, -1:dy:1 
grid = [ (x, y) for x in xs, y in ys ]
T_boundary = 200.
βs = [13221., 27., 5., 1.]
ω = 1.0

T₀ = T_boundary * ones(size(grid))
q₀ = map(xs -> compute_q(xs...), grid);

n = length(xs)
T_truth = deepcopy(T₀)

num_iters = 30
ε = zeros(num_iters);

## NEWTON ITERATIONS
#=================================================================#

for i in 1:num_iters
    R, ∂R∂T = compute_residual_and_grad(T_truth, βs, q₀, dx, dy)
    ΔT = newton_step(T_truth, R, ∂R∂T)
    ε[i] = maximum(abs.(ΔT))
    # display(R)
    println("Newton step error: $(ε[i])")
    T_truth .+= ω * ΔT
end

##
k_truth = compute_k.(βs..., T_truth)

##
using Plots
plotlyjs()

##
plot(ε)

##
surf1 = surface(xs, ys, q₀)

##
surf2 = surface(xs, ys, T_truth)

##
surf3 = surface(xs, ys, k_truth)

## ADJOINT ITERATIONS
#=================================================================#

function cost(T, T_truth, indices) 
    norm(T[index] - T_truth[index] for index in indices)
end

∇cost(T, T_truth, indices) = gradient(cost, T, T_truth, indices)


# indices = rand(1:length(T_truth[:]), (Int ∘ floor)(length(T_truth[:]) / 2))
indices = eachindex(T_truth)

##
β =  [13000., 10., 10., 5.]
γs = [100., 5., 5., 5.]

num_adjoint = 100
num_iters   = 30

T₀   = T_boundary * ones(size(grid))
q₀   = map(xs -> compute_q(xs...), grid);

T    = copy(T₀)
ε    = zeros(num_adjoint, num_iters)
εa   = zeros(num_adjoint)
β_l  = [β'; zeros(num_adjoint - 1, length(β)) ];
R    = zeros(size(T))
∂R∂T = zeros(size(R)..., size(T)...);

##
for i in 1:num_adjoint
    
    println("\nAdjoint Iteration — $i")
    display(β)

    for j in 1:num_iters
        R, ∂R∂T = compute_residual_and_grad(T, β, q₀, dx, dy)
        ΔT      = newton_step(T, R, ∂R∂T)
        ε[i,j]  = maximum(abs.(ΔT))
        # println("Newton step error: $(ε[i,j])")
        T     .+= ω * ΔT
    end

    εa[i] = cost(T, T_truth, indices)
    println("\nCost function value: $(εa[i])")

    # println("Computing df/dT")
    dfdT = ∇cost(T, T_truth, indices)[1]

    # println("Computing ∂f/∂x")
    ∂f∂x = zeros(length(β))    

    # println("Computing ∂R/∂x")
    ∂R∂x = ∇residual_β(T, β, q₀, dx, dy, T_boundary)

    # display(∂R∂x)

    # println("Solving adjoint problem: (∂R/∂T)ᵗ φ = -(df/dT)ᵗ")
    # φ = solve_adjoint(T, ∂R∂T, dfdT)
    # dfdβ = total_derivative_adjoint(∂f∂x, φ, ∂R∂x)

    dTdx = solve_direct(β, T, ∂R∂x, ∂R∂T)
    dfdβ = total_derivative_direct(∂f∂x, dTdx, dfdT)

    # display(dfdβ)

    # println("Updating design variables β")
    @. β -= γs * dfdβ
    β_l[i,:] = β' 
end


##
display(β)

##
k = compute_k.(β..., T)

##
using Plots
plotlyjs()

##
plot(εa, xlabel = "Iterations", ylabel = "Grid Error")

##
surface(xs, xs, T)
surface!(xs, xs, T_truth)

##
surface(xs, xs, k)
surface!(xs, xs, k_truth)

##
# residual_discrete_k(T, k, q, h) = first(k) * (Δ²x(T, h) + Δ²y(T, h)) + Δx(k, h) * Δx(T, h) + Δy(k, h) * Δy(T, h) + q
# ∇residual_discrete_k(Ts, βs, q, h) = gradient(residual_discrete_k, Ts, βs, q, h)

# function update_k_grid_ghost_points!(k_grid)
#     k_grid[1,:] = 2 * k_grid[2,:] .- k_grid[3,:]
#     k_grid[end,:] = 2 * k_grid[end-1,:] .- k_grid[end-2,:]
#     k_grid[:,1] = 2 * k_grid[:,2] .- k_grid[:,3]
#     k_grid[:,end] = 2 * k_grid[:,end-1] .- k_grid[:,end-2]

#     k_grid
# end
