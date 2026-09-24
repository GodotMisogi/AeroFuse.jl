module AeroFuseKernelAbstractionsExt

using AeroFuse
using KernelAbstractions
using LinearAlgebra
using StaticArrays

const KA = KernelAbstractions

import AeroFuse.PotentialFlow: AbstractPotentialFlowElement, device_solve_linear, device_induced_sum,
                               component_blocks, influence_coefficient, has_bc_override, apply_bc_row!,
                               boundary_condition, control_point, normal_vector

## Precision and transfer
#==========================================================================================#

"""
    device_float(backend)

Floating-point type used for device-side computation: `Float64` when the backend supports it,
otherwise `Float32` (e.g. Apple Metal).
"""
device_float(backend) = KA.supports_float64(backend) ? Float64 : Float32

# Rebuild a parametric element `E{S}` as `E{T}`, converting each floating field. Elements are
# plain structs of `SVector{3}`, scalars and `Bool` flags, so field-wise conversion is exact
# in structure and lets one routine cover every element type.
_convert_field(::Type{T}, x :: SVector{N}) where {T,N} = SVector{N,T}(x)
_convert_field(::Type{T}, x :: Bool) where T = x
_convert_field(::Type{T}, x :: Real) where T = T(x)

function convert_element(::Type{T}, el :: E) where {T, E <: AbstractPotentialFlowElement}
    wrapper = Base.typename(E).wrapper
    return wrapper{T}(ntuple(i -> _convert_field(T, getfield(el, i)), fieldcount(E))...)
end

# Upload a host collection to the backend with element-wise conversion.
function to_device(backend, host :: AbstractArray{E}) where E
    dev = KA.allocate(backend, E, length(host))
    copyto!(dev, vec(collect(host)))
    return dev
end

to_device(backend, ::Type{T}, elements :: AbstractArray{<: AbstractPotentialFlowElement}) where T =
    to_device(backend, map(el -> convert_element(T, el), vec(elements)))

to_device(backend, ::Type{T}, points :: AbstractArray{<: SVector{3}}) where T =
    to_device(backend, map(SVector{3,T}, vec(points)))

## Kernels
#==========================================================================================#

# Minimal collocation target so the device path reuses `influence_coefficient` verbatim,
# rather than restating the `velocity·n̂` formula.
struct Collocation{T} <: AbstractPotentialFlowElement
    rc     :: SVector{3,T}
    normal :: SVector{3,T}
end

# One column block of the AIC: every collocation row against one concrete-typed source block.
@kernel function aic_block_kernel!(A, @Const(targets), @Const(sources), col0)
    i, j = @index(Global, NTuple)
    @inbounds A[i, col0 + j] = influence_coefficient(sources[j], targets[i])
end

# Accumulate `f(r_i, el_j, s_j, V̂)` over one concrete-typed source block into `out[i]`.
# Each thread owns one output point, so accumulation is race-free across a launch.
@kernel function induced_sum_kernel!(out, f, @Const(points), @Const(sources), @Const(strengths), V_hat)
    i = @index(Global, Linear)
    @inbounds begin
        r = points[i]
        v = out[i]
        for j in eachindex(sources)
            v += f(r, sources[j], strengths[j], V_hat)
        end
        out[i] = v
    end
end

## Row-window adaptor for boundary-condition overrides
#==========================================================================================#

# Presents a host copy of rows `offset+1 : offset+size(rows,1)` of the full AIC under their
# global row indices, so element `apply_bc_row!` methods run unchanged on a partial download.
struct RowWindow{T, M <: AbstractMatrix{T}} <: AbstractMatrix{T}
    rows   :: M
    offset :: Int
    nrows  :: Int
end

Base.size(W :: RowWindow) = (W.nrows, size(W.rows, 2))
Base.getindex(W :: RowWindow, i :: Int, j :: Int) = W.rows[i - W.offset, j]
Base.setindex!(W :: RowWindow, v, i :: Int, j :: Int) = (W.rows[i - W.offset, j] = v)

function apply_bc_overrides!(AIC, boco, elements, blocks, U, Ups, Ω)
    N = length(elements)
    for (range, block) in blocks
        has_bc_override(first(block)) || continue
        rows = Array(AIC[range, :])
        W = RowWindow(rows, first(range) - 1, N)
        for i in range
            apply_bc_row!(W, boco, i, elements[i], elements, U, Ups, Ω)
        end
        AIC[range, :] .= to_matrix(KA.get_backend(AIC), rows)
    end
    return AIC, boco
end

function to_matrix(backend, host :: AbstractMatrix{T}) where T
    dev = KA.allocate(backend, T, size(host))
    copyto!(dev, host)
    return dev
end

## Backend implementations of the core hooks
#==========================================================================================#

function device_solve_linear(backend :: KA.Backend, elements, U, Ups, Ω)
    T = device_float(backend)
    N = length(elements)
    blocks = component_blocks(elements)

    # Collocation targets share one concrete type, so the rows need no per-block dispatch.
    targets = to_device(backend, map(el -> Collocation{T}(control_point(el), normal_vector(el)), vec(collect(elements))))

    AIC = KA.allocate(backend, T, N, N)
    kernel! = aic_block_kernel!(backend)
    for (range, block) in blocks
        sources = to_device(backend, T, block)
        kernel!(AIC, targets, sources, first(range) - 1; ndrange = (N, length(range)))
    end

    # Right-hand side assembled on the host: O(N), and it carries the prescribed onset fields.
    boco_host = boundary_condition(elements, U, Ups, Ω)
    apply_bc_overrides!(AIC, boco_host, elements, blocks, U, Ups, Ω)
    boco = to_device(backend, map(T, boco_host))

    # Strengths return to the host in the component layout of `boco_host`.
    strengths = similar(boco_host)
    copyto!(strengths, lu_solve(AIC, map(T, vec(collect(boco_host)))))

    return strengths, AIC, boco
end

# Factorize on the device (the O(N³) step), then apply the O(N²) triangular solves on the host.
# Device triangular solves are latency-bound (e.g. Metal's MPSMatrixSolveLU is ~50x slower than
# host BLAS at N ≈ 5000), whereas downloading the factors is a single cheap bulk copy.
function lu_solve(A, b :: Vector)
    F = lu(A; check = false)
    issuccess(F) || throw(SingularException(F.info))
    factors = F.factors isa Matrix ? F.factors : Array(F.factors)
    ipiv    = convert(Vector{LinearAlgebra.BlasInt}, Array(F.ipiv))
    return LU(factors, ipiv, F.info) \ b
end

function device_induced_sum(backend :: KA.Backend, f :: F, points, elements, strengths, V_hat) where F
    T = device_float(backend)
    pts = to_device(backend, T, points)
    out = KA.zeros(backend, SVector{3,T}, length(pts))
    kernel! = induced_sum_kernel!(backend)
    for (range, block) in component_blocks(elements)
        sources = to_device(backend, T, block)
        strs    = to_device(backend, map(T, vec(collect(strengths))[range]))
        kernel!(out, f, pts, sources, strs, SVector{3,T}(V_hat); ndrange = length(pts))
    end
    KA.synchronize(backend)
    return Array(out)
end

end
