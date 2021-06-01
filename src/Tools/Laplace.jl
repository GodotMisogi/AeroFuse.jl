module Laplace

using StaticArrays
using LinearAlgebra
using ..AeroMDAO: Point2D, Point3D

abstract type AbstractLaplace end

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

## 2D singularities
#============================================#

struct Singularity2D{T <: Real} <: AbstractLaplace
	str :: T
	r 	:: Point2D{T}
end

source_velocity(src :: Singularity2D, x, y) = SVector(src.str / (2π) * (x - src.r.y) / ((x - src.r.y)^2 + (y - src.r.x)^2), str / (2π) * (y - src.r.x) / ((x - src.r.y)^2 + (y - src.r.x)^2))
source_potential(src :: Singularity2D, x, y) = src.str / (4π) * log((x - src.r.y)^2 + (y - src.r.x)^2)
source_stream(src :: Singularity2D, x, y) = src.str / (2π) * atan(y - src.r.x, x - src.r.y)

doublet_velocity(dub :: Singularity2D, x, y) = SVector(dub.str / (2π) * ((x - dub.r.y)^2 - (y - dub.r.x)^2) / ((x - dub.r.y)^2 + (y - dub.r.x)^2)^2, - dub.str / (2π) * 2 * (x - dub.r.y) * (y - dub.r.x) / ((x - dub.r.y)^2 + (y - dub.r.x)^2)^2)
doublet_potential(dub :: Singularity2D, x, y) = -dub.str / (2π) * (y - dub.r.x) / ((x - dub.r.y)^2 + (y - dub.r.x)^2)
doublet_stream(dub :: Singularity2D, x, y) = -dub.str / (2π) * (y - dub.r.x) / ((x - dub.r.y)^2 + (y - dub.r.x)^2)

vortex_velocity(vor :: Singularity2D, x, y) = SVector(-vor.str / (2π) * (y - vor.r.x) / ((x - vor.r.y)^2 + (y - vor.r.x)^2), str / (2π) * (x - vor.r.y) / ((x - vor.r.y)^2 + (y - vor.r.x)^2))
vortex_potential(vor :: Singularity2D, x, y) = vor.str / (2π) * atan(y - vor.r.x, x - vor.r.y)
vortex_stream(vor :: Singularity2D, x, y) = -vor.str / (4π) * log((x - vor.r.y)^2 + (y - vor.r.x)^2)


struct Uniform2D{T <: Real} <: AbstractLaplace
	mag :: T
	ang :: T
	Uniform2D{T}(mag, ang) where T <: Real = new(mag, deg2rad(ang))
end

Uniform2D(mag :: T, ang :: T) where T <: Real = Uniform2D{T}(mag, ang)
Uniform2D(mag, ang) = Uniform2D(promote(mag, ang)...)

velocity(uni :: Uniform2D) = let (sa, ca) = sincos(uni.ang); uni.mag * SVector(ca, sa) end
potential(uni :: Uniform2D, x, y) = uni.mag * (x * cos(uni.ang) + y * sin(uni.ang))
stream(uni :: Uniform2D, x, y) = uni.mag * (y * cos(uni.ang) - x * sin(uni.ang))

## 3D singularities
#============================================#

# struct Singularity3D{T <: Real} <: AbstractLaplace
# 	str :: T
# 	r 	:: Point3D{T}
# end

# source_velocity(src :: Source2D, x, y, z)
# source_potential(src :: Source2D, x, y, z) 
# source_stream(src :: Source2D, x, y, z) 

## Freestream
#============================================#

struct Freestream{T <: Real} <: AbstractLaplace
	V :: T
	α :: T
	β :: T
	Ω :: SVector{3,T}
	Freestream{T}(V, α_deg, β_deg, Ω) where T <: Real = new(V, deg2rad(α_deg), deg2rad(β_deg), Ω)
end

"""
	Freestream(V, α, β, Ω)

A Freestream flow in spherical polar coordinates with magnitude ``V``, angle-of-attack ``α``, side-slip angle ``β``, and a quasi-steady rotation vector ``\\Omega``.
"""
Freestream(V, α_deg, β_deg, Ω :: AbstractVector{T}) where T <: Real = Freestream{T}(V, α_deg, β_deg, Ω)

"""
	Freestream(U, Ω)

A Freestream flow in Cartesian coordinates with vector ``U`` and quasi-steady rotation vector ``\\Omega``.
"""
Freestream(U :: AbstractVector{T}, Ω :: AbstractVector{T}) where T <: Real = Freestream{T}(cartesian_to_freestream(U)..., Ω)

"""
	freestream_to_cartesian(r, θ, φ)

Convert freestream flow (spherical polar) coordinates to Cartesian coordinates.
"""
freestream_to_cartesian(r, θ, φ) = r .* SVector(cos(θ) * cos(φ), -sin(φ), sin(θ) * cos(φ))

"""
	cartesian_to_freestream(U)

Convert Cartesian coordinates to freestream (spherical polar) flow coordinates.
"""
cartesian_to_freestream(U) = SVector(norm(U), -atand(U[3], U[1]), -atand(U[2] / √(U[1]^2 + U[3]^2)))

"""
	velocity(freestream :: Freestream)

Compute the velocity of a `Freestream`.
"""
velocity(freestream :: Freestream) = freestream_to_cartesian(freestream.V, freestream.α, freestream.β)

"""
	aircraft_velocity(freestream :: Freestream)

Compute the velocity of Freestream in the aircraft reference frame.
"""
aircraft_velocity(freestream :: Freestream) = -velocity(freestream)

end