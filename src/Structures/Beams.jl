module Beams

using SparseArrays
using StaticArrays

abstract type AbstractBeam end

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

# von_mises_stress(spar :: Tube) = 

# E * I / L^3 .* [ 12 .* [ 1. -1. ; -1.  1.]  6 .* [-1. -1.; 1. 1.] ; 
#                   6 .* [-1.  1. ; -1.  1.]  2 .* [ 2.  1.; 1. 2.] ]

## Coefficient matrices
J_coeffs(A, B, L) = A * B / L .* [ 1 -1; -1  1 ]

Iyy_coeffs(E, I, L) = let k0 = 12, k1 = 6L, k2 = 2L^2;
    E * I / L^3 .* @SMatrix [  k0 -k1 -k0 -k1 ;
                              -k1 2k2  k1  k2 ;
                              -k0  k1  k0  k1 ;
                              -k1  k2  k1 2k2 ] end

Izz_coeffs(E, I, L) = let k0 = 12, k1 = 6L, k2 = 2L^2;
    E * I / L^3 .* @SMatrix [  k0  k1 -k0  k1 ;
                               k1 2k2 -k1  k2 ;
                              -k0 -k1  k0 -k1 ;
                               k1  k2 -k1 2k2 ] end

## Composite deflection stiffness matrix of adjacent beam elements
function deflection_stiffness_matrix(Es, Is, Ls, direction = :z)
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
function torsional_stiffness_matrix(Es, Is, Ls)
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
function tube_stiffness_matrix(Es :: AbstractVector{T}, Gs :: AbstractVector{T}, As :: AbstractVector{T}, Iys :: AbstractVector{T}, Izs :: AbstractVector{T}, Js :: AbstractVector{T}, Ls :: AbstractVector{T}) where T <: Real
    # Check if sizes match
    @assert length(Es) == (length ∘ zip)(Gs, Iys, Izs, Js, Ls) "Lengths of coefficients must be the same as the dimension."

    Izs = deflection_stiffness_matrix(Es, Izs, Ls, :z)
    Iys = deflection_stiffness_matrix(Es, Iys, Ls, :y)
    As  = torsional_stiffness_matrix(Es, As, Ls)      
    Js  = torsional_stiffness_matrix(Gs, Js, Ls)

    K   = blockdiag(Iys, Izs, As, Js)
end

tube_stiffness_matrix(E :: Real, G :: Real, A :: Real, Iy :: Real, Iz :: Real, J :: Real, L :: Real, num :: Integer) = tube_stiffness_matrix(fill(E, num), fill(G, num), fill(A, num), fill(Iy, num), fill(Iz, num), fill(J, num), fill(L / num, num))

function tube_stiffness_matrix(x :: AbstractMatrix{<: Real}) 
    @assert size(x)[2] == 7 "Input must have 7 columns."
    @views tube_stiffness_matrix(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7])
end

function tube_stiffness_matrix(mat :: Material, tubes :: Vector{<: Tube})
    n    = length(tubes)
    As   = area.(tubes)
    Js   = polar_moment_of_inertia.(tubes)
    Iyys = moment_of_inertia.(tubes)
    Izzs = moment_of_inertia.(tubes)
    Ls   = length.(tubes)

    # Material propeties
    E = fill(elastic_modulus(mat), n)
    G = fill(shear_modulus(mat), n)

    tube_stiffness_matrix(E, G, As, Iyys, Izzs, Js, Ls)
end

end