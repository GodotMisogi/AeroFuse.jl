using LinearAlgebra
using Zygote
using StaticArrays
using Base.Iterators

## Finite differencing
central_2nd(Tm, T, Tp) = Tm - 2T + Tp
central_2nd_diff(Tm, T, Tp, h) = central_2nd(Tm, T, Tp) / h^2
central_1st_diff(Tm, Tp, h) = (Tm - Tp) / 2h

compute_k(β₁, β₂, T) = β₁/T + β₂

∇compute_k(β₁, β₂, T) = gradient(compute_k, β₁, β₂, T)

compute_q(x, y) = 5e5 * (1 - (x - 0.5)^2) * (1 - (y - 0.5)^2)

function get_T(T_grid, i, j, T_boundary)
    if i < 1 || j < 1 || i > length(T_grid[:,1])|| j > length(T_grid[:,1])
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

function residual(T, βs, q, h)
    k = compute_k(βs..., first(T)) 
    ∂k∂T = ∇compute_k(βs..., first(T))[3]
    
    k * (Δ²x(T, h) + Δ²y(T, h)) + ∂k∂T * (Δx(T, h)^2 + Δy(T, h)^2 ) + q
end

# function residual_boundary(T, T_boundary, βs, q, h)
#     k = compute_k(βs..., first(T)) 
#     ∂k∂T = ∇compute_k(βs..., first(T))[3]
    
# end

∇residual(Ts, βs, q, h) = gradient(residual, Ts, βs, q, h)

## Grid setup
Δh = 0.2
xs = 0:Δh:1
grid = product(xs, xs)
T_boundary = 200.
βs = [13221, 27]
ω = 0.8

T₀ = T_boundary * ones(size(grid))
q₀ = map(xs -> compute_q(xs...), grid);

n = length(xs)
T = deepcopy(T₀)
R = zeros(size(T)...)
∂R∂T = zeros(n, n, n, n)
num_iters = 15
ε = zeros(num_iters)

## Iteration
for i in 1:num_iters
    for R_ind in CartesianIndices(R)
        k, l = R_ind.I
        Ts = stencil_1(T, k, l, T_boundary)
        
        # Residual 
        R[R_ind] = residual(Ts, βs, q₀[R_ind], Δh)
        
        # Jacobian
        grad = ∇residual(Ts, βs, q₀[R_ind], Δh)[1]

        ∂R∂T[R_ind,R_ind] = grad[1]
        if k > 1 ∂R∂T[R_ind,k-1,l] = grad[2] end
        if k < n ∂R∂T[R_ind,k+1,l] = grad[3] end
        if l > 1 ∂R∂T[R_ind,k,l-1] = grad[4] end
        if l < n ∂R∂T[R_ind,k,l+1] = grad[5] end
    end
    
    # display(T)
    # display(R)
    # display(reshape(∂R∂T, size(T).^2))

    # Solve system
    ΔT = - reshape(reshape(∂R∂T, size(T).^2) \ R[:], size(T)) 
    ε[i] = maximum(abs.(ΔT))
    T .+= ω * ΔT
end

k = compute_k.(βs..., T)

##
using Plots
plotlyjs()

##
plot(ε)

##
surf1 = surface(xs, xs, q₀)

##
surf2 = surface(xs, xs, T)

##
surf3 = surface(xs, xs, k)

# begin
#     ΔT = ∂R∂T \ -ΔR
#     T₀ += ω * reshape(ΔT, size(T₀))
#     residual = max(abs.(ΔT))
# end#


# residual_discrete_k(T, k, q, h) = 
#     first(k) * (Δ²x(T) + Δ²y(T)) / h^2 + (Δx(k) / h) * (Δx(T) / h) + (Δy(k) / h) * (Δy(T) / h) + q

# ∇residual_discrete_k(Ts, βs, q, h) = gradient(residual_discrete_k, Ts, βs, q, h)

# function update_k_grid_ghost_points!(k_grid)
#     k_grid[1,:] = 2 * k_grid[2,:] .- k_grid[3,:]
#     k_grid[end,:] = 2 * k_grid[end-1,:] .- k_grid[end-2,:]
#     k_grid[:,1] = 2 * k_grid[:,2] .- k_grid[:,3]
#     k_grid[:,end] = 2 * k_grid[:,end-1] .- k_grid[:,end-2]

#     k_grid
# end
