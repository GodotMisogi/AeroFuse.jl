module Beams

## Package imports
#==========================================================================================#

using SparseArrays
using StaticArrays
using LinearAlgebra

abstract type AbstractBeam end

## Material definition
#==========================================================================================#

struct Material{T}
    elastic_modulus :: T
    shear_modulus   :: T
    yield_stress    :: T
    density         :: T
end

"""
    Material(E, J, σ_max, ρ)

Define a `Material` with positive, real elastic modulus ``E``, shear modulus ``G``, yield stress ``\\sigma_\\max``, and density ``\\rho``.
"""
function Material(E :: T, G :: T, σ_max :: T, ρ :: T) where T <: Real
    @assert E > 0. && G > 0. && σ_max > 0. && ρ > 0. "All quantities must be positive."
    Material{T}(E, G, σ_max, ρ)
end

elastic_modulus(mat :: Material) = mat.elastic_modulus
shear_modulus(mat :: Material)   = mat.shear_modulus
yield_stress(mat :: Material)    = mat.yield_stress
density(mat :: Material)         = mat.density


## Tube definition
#==========================================================================================#

struct Tube{T <: Real} <: AbstractBeam
    material  :: Material{T}
    length    :: T
    radius    :: T
    thickness :: T
    function Tube(material :: Material{T}, length :: T, radius :: T, thickness :: T) where T <: Real
        @assert length > 0. && radius > 0. && thickness > 0. "Length, outer radius and thickness must be positive."
        new{T}(material, length, radius, thickness)
    end
end

Base.length(tube :: Tube) = tube.length
material(tube :: Tube)    = tube.material
radius(tube :: Tube)      = tube.radius
thickness(tube :: Tube)   = tube.thickness

radii(tube :: Tube) = radius(tube) - thickness(tube), radius(tube)

function area(tube :: Tube)
    r1, r2 = radii(tube)
    π * (r2^2 - r1^2)
end

function moment_of_inertia(tube :: Tube)
    r1, r2 = radii(tube)
    π/4 * (r2^4 - r1^4)
end

function polar_moment_of_inertia(tube :: Tube)
    r1, r2 = radii(tube)
    π/2 * (r2^4 - r1^4)
end

## Stress calculations (NEEDS CHECKING)
#==========================================================================================#

principal_stress(E, L, R, dx, dθ_yz) = E * (dx / L + R * dθ_yz / L)
torsional_stress(G, L, R, dθ_x)    = G * R * dθ_x / L
von_mises_stress(σ_xx, σ_xt)      = √(σ_xx^2 + 3σ_xt^2)

function von_mises_stress(tube :: Tube, ds, θs)
    E = (elastic_modulus ∘ material)(tube)
    G = (shear_modulus ∘ material)(tube)
    R = radius(tube)
    L = length(tube)

    dx     = ds[1]
    dθ_yz  = norm(θs[2:end])
    σ_xx_1 = principal_stress(E, L, R,  dx, dθ_yz)
    σ_xx_2 = principal_stress(E, L, R, -dx, dθ_yz)
    σ_xt   = torsional_stress(G, L, R, θs[1])

    von_mises_stress.(SVector(σ_xx_1, σ_xx_2), σ_xt)
end

## Stiffness matrices
#==========================================================================================#

# E * I / L^3 .* [ 12 .* [ 1. -1. ; -1.  1.]  6 .* [-1. -1.; 1. 1.] ;
#                   6 .* [-1.  1. ; -1.  1.]  2 .* [ 2.  1.; 1. 2.] ]

## Coefficient matrices
J_coeffs(A, B, L) = A * B / L .* [ 1 -1; -1  1 ]

k_coeffs(k0, k1, k2) = @SMatrix [  k0 -k1 -k0 -k1 ;
                                  -k1 2k2  k1  k2 ;
                                  -k0  k1  k0  k1 ;
                                  -k1  k2  k1 2k2 ]

Iyy_coeffs(E, I, L) = let k0 = 12, k1 =  6L, k2 = 2L^2; E * I / L^3 * k_coeffs(k0, k1, k2) end
Izz_coeffs(E, I, L) = let k0 = 12, k1 = -6L, k2 = 2L^2; E * I / L^3 * k_coeffs(k0, k1, k2) end

## Composite deflection stiffness matrix of adjacent beam elements
function bending_stiffness_matrix(Es, Is, Ls, direction = :z)
    n = length(Es)
    mat = @MMatrix zeros(2(n+1), 2(n+1))
    block_size = CartesianIndices((4,4))

    for (c_ind, E, I, L) in zip(CartesianIndex.(0:2:2(n-1), 0:2:2(n-1)), Es, Is, Ls)
        coeffs = ifelse(direction == :z, Izz_coeffs(E, I, L), Iyy_coeffs(E, I, L))
        @. mat[c_ind + block_size] += coeffs
    end

    sparse(mat)
end

## Composite torsional stiffness matrix for adjacent beam elements
function axial_stiffness_matrix(Es, Is, Ls)
    n = length(Es)
    mat = @MMatrix zeros(n+1, n+1)
    block_size = CartesianIndices((2,2))

    for (c_ind, E, I, L) in zip(CartesianIndex.(0:n-1, 0:n-1), Es, Is, Ls)
        coeffs = J_coeffs(E, I, L)
        @. mat[c_ind + block_size] += coeffs
    end

    sparse(mat)
end

## Full stiffness matrix for one element
tube_stiffness_matrix(E :: Real, G :: Real, A :: Real, Iy :: Real, Iz :: Real, J :: Real, L :: Real) = blockdiag((sparse ∘ Iyy_coeffs)(E, Iy, L), (sparse ∘ Izz_coeffs)(E, Iz, L), (sparse ∘ J_coeffs)(E, A, L), (sparse ∘ J_coeffs)(G, J, L))

## Composite stiffness matrix for adjacent tube finite elements
function tube_stiffness_matrix(Es :: AbstractVector, Gs :: AbstractVector, As :: AbstractVector, Iys :: AbstractVector, Izs :: AbstractVector, Js :: AbstractVector, Ls :: AbstractVector)
    # Check if sizes match
    @assert length(Es) == (length ∘ zip)(Gs, Iys, Izs, Js, Ls) "Lengths of coefficients must be the same as the dimension."

    Izs = bending_stiffness_matrix(Es, Izs, Ls, :z)
    Iys = bending_stiffness_matrix(Es, Iys, Ls, :y)
    As  = axial_stiffness_matrix(Es, As, Ls)
    Js  = axial_stiffness_matrix(Gs, Js, Ls)

    blockdiag(Iys, Iys, As, Js)
end

tube_stiffness_matrix(E :: Real, G :: Real, A :: Real, Iy :: Real, Iz :: Real, J :: Real, L :: Real, num :: Integer) = tube_stiffness_matrix(fill(E, num), fill(G, num), fill(A, num), fill(Iy, num), fill(Iz, num), fill(J, num), fill(L / num, num))

tube_stiffness_matrix(tube :: Tube) = tube_stiffness_matrix(elastic_modulus(material(tube)), shear_modulus(material(tube)), area(tube), moment_of_inertia(tube), moment_of_inertia(tube), polar_moment_of_inertia(tube), length(tube))

function tube_stiffness_matrix(x :: AbstractMatrix)
    @assert size(x)[2] == 7 "Input must have 7 columns."
    @views tube_stiffness_matrix(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7])
end

tube_stiffness_matrix(tubes :: Vector{<: Tube}) = tube_stiffness_matrix((elastic_modulus ∘ material).(tubes), (shear_modulus ∘ material).(tubes), area.(tubes), moment_of_inertia.(tubes), moment_of_inertia.(tubes), polar_moment_of_inertia.(tubes), length.(tubes))

function build_stiffness_matrix(Ks, constraint_indices)
    # Temporary reshaping for sparse matrix construction
    num_Ks = length(Ks)
    D      = @views reshape(reduce(hcat, Ks), 12, 12, num_Ks)

    # First element
    D_start = @views D[1:6,1:6,1]

    # Last element
    D_end   = @views D[7:end,7:end,end]

    # Summing adjacent elements corresponding to diagonals
    D_mid   = @views D[7:end,7:end,1:end-1] + D[1:6,1:6,2:end]

    # Upper and lower diagonals
    D_12    = @views D[1:6,7:end,:]
    D_21    = @views D[7:end,1:6,:]

    # Build sparse matrix
    stiffness = spzeros(6 * (num_Ks + 2), 6 * (num_Ks + 2))
    @views stiffness[7:12,7:12]           = D_start
    @views stiffness[end-5:end,end-5:end] = D_end

    for m in 1:num_Ks
        mid_inds = 6(m+1)+1:6(m+1)+6
        off_inds = 6(m)+1:6(m)+6
        if m < num_Ks
            @views stiffness[mid_inds,mid_inds] = D_mid[:,:,m]
        end
        @views stiffness[off_inds,mid_inds] = D_12[:,:,m]
        @views stiffness[mid_inds,off_inds] = D_21[:,:,m]
    end

    # Fixed boundary condition by constraining the locations.
    cons = 6 * (constraint_indices .+ 1)
    arr = 1:6

    col = arr
    row = reduce(hcat, con .+ arr for con in cons)

    @views stiffness[CartesianIndex.(col, row)] .= 1e9
    @views stiffness[CartesianIndex.(row, col)] .= 1e9

    stiffness
end

## Cantilever setup
function solve_cantilever_beam(Ks, loads, constraint_indices)
    # Create the stiffness matrix from the array of individual stiffnesses
    # Also specifies the constraint location
    K = build_stiffness_matrix(Ks, constraint_indices)

    # Build force vector with constraint
    f = [ zeros(6); loads[:] ]

    # Solve FEM system
    x = K \ f

    # Throw away the junk values for the constraint
    reshape(x[7:end], 6, length(Ks) + 1)
end

end