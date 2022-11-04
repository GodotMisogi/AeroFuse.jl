module Laplace

## Package imports
#==========================================================================================#

using StaticArrays
using LinearAlgebra
using CoordinateTransformations, Rotations
using Setfield

import Base: *

const cfs = CartesianFromSpherical()
const sfc = SphericalFromCartesian()

import ..MathTools: Point2D, Point3D, magnitude, angle

import ..AeroFuse: velocity

abstract type AbstractLaplace end

## Legacy (to be removed?)
#==========================================================================================#

# Performs velocity and potential computations for an object on a grid
grid_data(object :: AbstractLaplace, xs) = velocity(object, xs), potential(object, xs)
velocity(object :: AbstractLaplace, xs) = map(x -> velocity(object, x...), xs) 
potential(object :: AbstractLaplace, xs) = map(x -> potential(object, x...), xs)

# Performs velocity and potential calculations on a grid
function grid_data(objects :: Vector{<: AbstractLaplace}, xs)
    vels = foldl((v1, v2) -> [ u .+ v for (u, v) ∈ zip(v1, v2) ], velocity(object, xs) for object ∈ objects)
    pots = foldl((v1, v2) -> v1 + v2, potential(object, xs) for object ∈ objects)
    
    vels, pots
end

# Generic attempt
abstract type AbstractSingularity <: AbstractLaplace end

# Generic geometric type
struct GeometricSingularity{T <: Number, K <: AbstractSingularity, P <: Integer, M <: Integer, N <: Integer} <: AbstractLaplace
    type :: K # Type of singularity (doublet, source, vortex, uniform, etc. if any more exist)
    strs :: Union{T,SVector{P,T}} # Represents P (constant = 1, linear = 2, quadratic = 3, etc.) singularity strengths
    rs   :: SMatrix{M,N,T} # Represents M points in N dimensions
end

## 2D singularities
#============================================#

abstract type AbstractSingularity2D <: AbstractSingularity end

# Geometric type
struct PointSingularity2D{T <: Number, N <: AbstractSingularity2D} <: AbstractLaplace
    type     :: N
    strength :: T
    x        :: T
    y        :: T
end

# Promoted-type constructor
function PointSingularity2D(type :: N, str, x, y) where N <: AbstractSingularity2D
    T = promote_type(eltype(str), eltype(x), eltype(y))
    PointSingularity2D{T,N}(type, str, x, y)
end

# Getters
strength(s :: PointSingularity2D) = s.strength
x(s :: PointSingularity2D) = s.x
y(s :: PointSingularity2D) = s.y

# Dispatchers
potential(src :: PointSingularity2D, x, y) = potential(src.type, src.strength, x - src.x, y - src.y) 
velocity(src :: PointSingularity2D, x, y) = velocity(src.type, src.strength, x - src.x, y - src.y) 
stream(src :: PointSingularity2D, x, y) = stream(src.type, src.strength, x - src.x, y - src.y) 

# Scaling law
# Base.*(a :: Real, src :: PointSingularity2D) = @set src.strength = src.strength * a 

# Linearity by generated functions?
# @generated velocity()

# Source relations
struct Source2D <: AbstractSingularity2D end

Source2D(str, x, y) = PointSingularity2D(Source2D(), str, x, y)
Source2D(str, r) = PointSingularity2D(Source2D(), str, r[1], r[2])

potential(:: Source2D, σ, x, y) = σ / 4π * log(x^2 + y^2)
velocity(:: Source2D, σ, x, y) = @SVector [ σ / (2π) * x / (x^2 + y^2), str / (2π) * y / (x^2 + y^2) ]
stream(:: Source2D, σ, x, y) = σ / 2π * atan(y, x)

# Doublet relations
struct Doublet2D <: AbstractSingularity2D end

Doublet2D(str, x, y) = PointSingularity2D(Doublet2D(), str, x, y)
Doublet2D(str, r) = PointSingularity2D(Doublet2D(), str, r[1], r[2])

potential(:: Doublet2D, μ, x, y) = -μ / (2π) * y / (x^2 + y^2)
velocity(:: Doublet2D, μ, x, y) = @SVector [μ / (2π) * (x^2 - y^2) / (x^2 + y^2)^2, - μ / (2π) * 2 * x * y / (x^2 + y^2)^2]
stream(:: Doublet2D, μ, x, y) = -μ / (2π) * y / (x^2 + y^2)

# Vortex relations
struct Vortex2D <: AbstractSingularity2D end

Vortex2D(str, x, y) = PointSingularity2D(Vortex2D(), str, x, y)
Vortex2D(str, r) = PointSingularity2D(Vortex2D(), str, r[1], r[2])

potential(:: Vortex2D, γ, x, y) = γ / (2π) * atan(y, x)
velocity(:: Vortex2D, γ, x, y) = @SVector [-γ / (2π) * y / (x^2 + y^2), str / (2π) * x / (x^2 + y^2) ]
stream(:: Vortex2D, γ, x, y) = -γ / (4π) * log(x^2 + yp^2)

# Uniform relations

# struct Uniform2D <: AbstractSingularity2D end

# Uniform2D(str, r) = PointSingularity2D(Uniform2D(), str, r)

struct Uniform2D{T <: Real} <: AbstractSingularity2D
    magnitude :: T
    angle     :: T
    Uniform2D{T}(mag, ang) where T <: Real = new(mag, deg2rad(ang))
end

magnitude(uni :: Uniform2D) = uni.magnitude
angle(uni :: Uniform2D) = uni.angle

Uniform2D(mag :: T, ang :: T) where T <: Real = Uniform2D{T}(mag, ang)
Uniform2D(mag, ang) = Uniform2D(promote(mag, ang)...)
Uniform2D(velocity) = Uniform2D((sqrt ∘ sum)(velocity.^2), atan(velocity[2], velocity[1]))
Uniform2D(; angle)  = Uniform2D(1., angle)

potential(uni :: Uniform2D, x, y) = let (sa, ca) = sincos(uni.angle); uni.magnitude * (x * ca + y * sa) end
velocity(uni :: Uniform2D) = let (sa, ca) = sincos(uni.angle); uni.magnitude * SVector(ca, sa) end
stream(uni :: Uniform2D, x, y) = let (sa, ca) = sincos(uni.angle); uni.magnitude * (y * ca - x * sa) end

## 3D singularities
#============================================#

# Traits
abstract type AbstractSingularity3D <: AbstractSingularity end

# Geometric type
struct PointSingularity3D{T <: Number, N <: AbstractSingularity3D} <: AbstractLaplace
    type     :: N
    strength :: T
    r        :: SVector{3,T}
end

# Promoted-type constructor
function PointSingularity3D(type :: N, str, r) where N <: AbstractSingularity3D
    T = promote_type(eltype(str), eltype(r))
    PointSingularity3D{T,N}(type, str, r)
end

# Getters
strength(s :: PointSingularity3D) = s.strength
position(s :: PointSingularity3D) = s.r

# Dispatchers (streamfunction is not uniquely defined in 3D)
potential(src :: PointSingularity3D, r) = potential(src.type, src.strength, r - src.r) 
velocity(src :: PointSingularity3D, r) = velocity(src.type, src.strength, r - src.r) 

# Point source
struct Source3D <: AbstractSingularity3D end

Source3D(str, x, y, z) = PointSingularity3D(Source3D(), str, @SVector [x, y, z])
Source3D(str, r) = PointSingularity3D(Source3D(), str, r)

potential(:: Source3D, σ, r) = σ / (4π * norm(r))
velocity(:: Source3D, σ, r) = σ / 4π * r / norm(r)^3

# Doublet/dipole aligned along arbitrary axis
struct Doublet3D{T <: Number} <: AbstractSingularity3D 
    p :: SVector{3,T} # Orientation axis of the dipole
end

Doublet3D(str, x, y, z, p = SVector(1., 0., 0.)) = PointSingularity3D(Doublet3D(normalize(p)), str, @SVector [x, y, z])
Doublet3D(str, r, p = SVector(1., 0., 0.)) = PointSingularity3D(Doublet3D(normalize(p)), str, r)

function velocity(src :: Doublet3D, μ, r)
    x, y, z = r
    rn = norm(r)
    cφ = dot(src.p, r)

    μ / (4π * rn^3) * @SVector [ 3z * x / rn^2, 3z * y / rn^2, 3 * (cφ^2 - 1) ]
end

potential(dub :: Doublet3D, μ, r) = -μ / 4π * dot(dub.p, r) / norm(r)^3

struct ConstantStrengthLineSingularity3D{T <: Number, K <: AbstractSingularity3D} <: AbstractSingularity
    type     :: K
    strength :: T
    r1       :: SVector{3,T}
    r2       :: SVector{3,T}
end

# Alias
const CSLine3D = ConstantStrengthLineSingularity3D

(l :: CSLine3D)(; type = l.type, strength = l.strength, r1 = l.r1, r2 = l.r2) = CSLine3D(type, strength, r1, r2)

# function Base.show(io :: IO, src :: ConstantStrengthLine3D)
#     println(io, "Constant strength 3D potential line")
#     println(io, "    (r1, r2) = ", src.r1, " ", src.r2)
# end

function ConstantStrengthLineSingularity3D(type :: K, str, r1, r2, η) where K <: AbstractSingularity3D
    T = promote_type(typeof(str), eltype(r1), eltype(r2), eltype(η))
    CSLine{T,K}(type, str, r1, r2, η)
end

# Dispatchers (streamfunction is not uniquely defined in 3D)
potential(src :: CSLine3D, r) = potential(src.type, src.strength, r - src.r2) - potential(src.type, src.strength, r - src.r1)
velocity(src :: CSLine3D, r) = velocity(src.type, src.strength, r - src.r2) - velocity(src.type, src.strength, r - src.r1)

# Source line
SourceLine3D(str, r1, r2) = CSLine3D(Source3D(), str, r1, r2)

# Doublet line
DoubletLine3D(str, r1, r2, p = @SVector [1., 0., 0.]) = CSLine3D(Doublet3D(p), str, r1, r2)


## Freestream
#============================================#

abstract type AbstractFreestream end

"""
    Freestream(α, β, Ω)
    Freestream(α, β, Ω_x, Ω_y, Ω_z)
    Freestream(U, Ω)
    Freestream(; 
        alpha = 0., 
        beta  = 0., 
        omega = [0,0,0]
    )

Define freestream conditions with angle of attack ``α`` (degrees), sideslip angle ``β`` (degrees), and (quasi-steady) rotation vector ``Ω`` for a vortex lattice analysis.

Alternatively, provide the velocity vector ``U``, which is normalized to determine the angles.
"""
struct Freestream{T} <: AbstractFreestream
    alpha :: T
    beta  :: T
    omega :: SVector{3,T}
end

Base.broadcastable(fs :: Freestream) = Ref(fs)

Freestream(α_deg, β_deg, Ω) = let T = promote_type(eltype(α_deg), eltype(β_deg), eltype(Ω)); Freestream{T}(deg2rad(α_deg), deg2rad(β_deg), Ω) end

Freestream(α_deg, β_deg, Ω_x, Ω_y, Ω_z) = Freestream(deg2rad(α_deg), deg2rad(β_deg), SVector(Ω_x, Ω_y, Ω_z))

Freestream(U, Ω) = let (_, α, β) = cartesian_to_freestream(normalize(U)); Freestream(α, β, Ω) end

Freestream(; alpha = 0., beta = 0., omega = zeros(3)) = Freestream(alpha, beta, omega)

function Base.show(io :: IO, fs :: Freestream)
    println(io, "Freestream: ")
    for fname in fieldnames(typeof(fs))
        println(io, "    ", fname, " = ", getfield(fs, fname))
    end
end

potential(fs :: Freestream, r) = dot(velocity(fs), r)

"""
    velocity(freestream :: Freestream)

Compute the velocity of a `Freestream`.
"""
velocity(fs :: Freestream) = freestream_to_cartesian(1., fs.alpha, fs.beta)

"""
    freestream_to_cartesian(r, θ, φ)

Convert freestream flow (spherical polar) coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r * SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
    cartesian_to_freestream(U)

Convert Cartesian coordinates to freestream (spherical polar) flow coordinates.
"""
cartesian_to_freestream(U) = SVector(norm(U), -atand(U[3], U[1]), -atand(U[2], √(U[1]^2 + U[3]^2)))

# 2D versions
cartesian_to_freestream(u, w) = magnitude(u, w), angle(u, w)
freestream_to_cartesian(V, α) = V * cos(α), V * sin(α)

end