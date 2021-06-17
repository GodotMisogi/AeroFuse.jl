module Beams

using SparseArrays

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
function Material(E, G, σ_max, ρ)
    @assert E > 0. && G > 0. && σ_max > 0. && ρ > 0. "All quantities must be positive."
    Material{T}(E, G, σ_max, ρ)
end

struct Beam{T} <: AbstractBeam
    area   :: T
    Iy     :: T
    Iz     :: T
    length :: T
end

function Beam(A :: T, Iy :: T, Iz :: T, L :: T) where T <: Real
    @assert A > 0. && Iy > 0. && Iz > 0. && L > 0. "All quantities must be positive."
    Beam{T}(A, Iy, Iz, L)
end

struct Tube{T} <: AbstractBeam
    beam      :: Beam{T}
    J         :: T
    radius    :: T
    thickness :: T
end

function Tube(beam :: Beam{T}, radius :: T, thickness :: T) where T <: Real 
    @assert radius > 0. && thickness > 0. "Radius and thickness must be positive."   
    Tube{T}(beam, radius, thickness)
end

radii(tube :: Tube) = tube.radius - tube.thickness, tube.radius

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
J_coeffs(A, B, L) = A * B / L .* [  1 -1 ;
                                   -1  1 ]

Iy_coeffs(E, I, L) = let k0 = 12, k1 = 6L, k2 = 2L^2;
    E * I / L^3 .* [  k0 -k1 -k0 -k1 ;
                     -k1 2k2  k1  k2 ;
                     -k0  k1  k0  k1 ;
                     -k1  k2  k1 2k2 ] end

Iz_coeffs(E, I, L) = let k0 = 12, k1 = 6L, k2 = 2L^2;
    E * I / L^3 .* [  k0  k1 -k0  k1 ;
                      k1 2k2 -k1  k2 ;
                     -k0 -k1  k0 -k1 ;
                      k1  k2 -k1 2k2 ] end

## Composite deflection stiffness matrix of adjacent beam elements
function deflection_stiffness_matrix(Es, Is, Ls, direction = :z)
    n = length(Es)
    mat = zeros(2(n+1), 2(n+1))
    block_size = CartesianIndices((4,4))

    for (c_ind, E, I, L) in zip(CartesianIndex.(0:2:2(n-1), 0:2:2(n-1)), Es, Is, Ls)
        coeffs = ifelse(direction == :z, Iz_coeffs(E, I, L), Iy_coeffs(E, I, L))
        @. mat[c_ind + block_size] += coeffs
    end

    sparse(mat)
end

## Composite torsional stiffness matrix for adjacent beam elements
function torsional_stiffness_matrix(Es, Is, Ls)
    n = length(Es)
    mat = zeros(n+1, n+1)
    block_size = CartesianIndices((2,2))

    for (c_ind, E, I, L) in zip(CartesianIndex.(0:n-1, 0:n-1), Es, Is, Ls)
        coeffs = J_coeffs(E, I, L)
        @. mat[c_ind + block_size] += coeffs
    end

    sparse(mat)
end

## Full stiffness matrix for one element
tube_stiffness_matrix(E :: Real, G :: Real, A :: Real, Iy :: Real, Iz :: Real, J :: Real, L :: Real) = blockdiag((sparse ∘ Iy_coeffs)(E, Iy, L), (sparse ∘ Iz_coeffs)(E, Iz, L), (sparse ∘ J_coeffs)(E, A, L), (sparse ∘ J_coeffs)(G, J, L))

## Composite stiffness matrix for adjacent tube finite elements
function tube_stiffness_matrix(Es :: Vector{T}, Gs :: Vector{T}, As :: Vector{T}, Iys :: Vector{T}, Izs :: Vector{T}, Js :: Vector{T}, Ls :: Vector{T}) where T <: Real
    # Check if sizes match
    @assert length(Es) == (length ∘ zip)(Gs, Iys, Izs, Js, Ls) "Lengths of coefficients must be the same as the dimension."

    Izs = deflection_stiffness_matrix(Es, Izs, Ls, :z)
    Iys = deflection_stiffness_matrix(Es, Iys, Ls, :y)
    As  = torsional_stiffness_matrix(Es, As, Ls)      
    Js  = torsional_stiffness_matrix(Gs, Js, Ls)

    K   = blockdiag(Iys, Izs, As, Js)
end


end