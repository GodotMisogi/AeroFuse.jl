using LinearAlgebra
using Zygote
using StaticArrays
using Base.Iterators

##

# Finite differencing
central_2nd(T) = T[1] - 2 * T[2] + T[3]
central_2nd_diff(T, h) = central_2nd(T) / h^2
# central_2nd_diff(T, h) = diff(diff(T, dims = 2), dims = 2) / h^2
central_1st_diff(T, h) = (T[3] - T[1]) / (2h)^2

# [β₁, β₂, T]
compute_k(β₁, β₂, T) = β₁/T + β₂

∇compute_k(β₁, β₂, T) = gradient(compute_k, β₁, β₂, T)

compute_q(x, y) = 5e5 * (1 - (x - 0.5)^2) * (1 - (y - 0.5)^2)

function residual(Ts, βs, q, h)
    T_row = SVector(Ts[2], Ts[1], Ts[3]) # T_i-1, j   T_i,j T_i+1,j
    T_col = SVector(Ts[4], Ts[1], Ts[5])
    k = compute_k(βs..., Ts[1]) 
    ∂k∂T = ∇compute_k(βs..., Ts[1])[3]
    
    k * (central_2nd_diff(T_row, h) + central_2nd_diff(T_col, h)) + ∂k∂T * (central_1st_diff(T_row, h)^2 + central_1st_diff(T_col, h)^2) + q
end

∇residual(Ts, βs, q, h) = gradient(residual, Ts, βs, q, h)

function residual_discrete_k(Ts, ks, q, h)
    T_row = SVector(Ts[2], Ts[1], Ts[3]) # T_i-1, j   T_i,j T_i+1,j
    T_col = SVector(Ts[4], Ts[1], Ts[5])

    k_row = SVector(ks[2], ks[1], ks[3]) # k_i-1, j   k_i,j k_i+1,j
    k_col = SVector(ks[4], ks[1], ks[5])

    k_row[1] * (central_2nd_diff(T_row, h) + central_2nd_diff(T_col, h)) +
    central_1st_diff(k_row, h) * central_1st_diff(T_row, h) + central_1st_diff(k_col, h) * central_1st_diff(T_col, h) + q
end

∇residual_discrete_k(Ts, βs, q, h) = gradient(residual_discrete_k, Ts, βs, q, h)

function update_k_grid_ghost_points!(k_grid)
    k_grid[1, :] = 2 * k_grid[2, :] .- k_grid[3, :]
    k_grid[end, :] = 2 * k_grid[end-1, :] .- k_grid[end-2, :]
    k_grid[:, 1] = 2 * k_grid[:, 2] .- k_grid[:, 3]
    k_grid[:, end] = 2 * k_grid[:, end-1] .- k_grid[:, end-2]

    k_grid
end

@Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
@Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)


## Grid setup
Δh = 0.05
xs = 0:Δh:1
grid = product(xs, xs)
T_boundary = 200.
βs = [ 0.1, 0.5 ]

T₀ = T_boundary * ones(size(grid))
q₀ = map(xs -> compute_q(xs...), grid)

stencil_1(T, i, j) = SVector(T[i, j], T[i-1, j], T[i+1, j], T[i, j-1], T[i, j+1])

R = [  
        let Ts = stencil_1(T₀, indices...);
        residual(Ts, βs, q₀[indices...], Δh) end
        for indices in product(2:length(xs)-1, 2:length(xs)-1) 
    ]

∇R = 
    [  
        let Ts = stencil_1(T₀, indices...);
        ∇residual(Ts, βs, q₀[indices...], Δh)[1] end
        for indices in product(2:length(xs)-1, 2:length(xs)-1) 
    ]

# ∂R∂T = ???

# begin
#     ΔT = ∂R∂T \ -ΔR
#     T₀ += ω * reshape(ΔT, size(T₀))
#     residual = max(abs.(ΔT))
# end