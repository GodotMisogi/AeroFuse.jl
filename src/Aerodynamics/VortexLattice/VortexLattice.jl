module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations, CoordinateTransformations
using SplitApplyCombine
using TimerOutputs
using LabelledArrays
using Accessors
using PrettyTables
using ForwardDiff: jacobian!
using DiffResults: JacobianResult, jacobian, value

## Package imports
#==========================================================================================#

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Some tools
import ..Laplace: AbstractFreestream, Freestream

import ..AeroFuse: velocity, solve_linear, solve_nonlinear, solve_nonlinear!, surface_velocities, surface_coefficients

## Vortex types and methods
#==========================================================================================#

include("vortices.jl")

## Reference frames
#==========================================================================================#

## Axis transformations with respect to freestream
abstract type AbstractAxisSystem end

struct Geometry <: AbstractAxisSystem end
struct Body <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind <: AbstractAxisSystem end

"""
    velocity(freestream :: Freestream, ::Geometry)

Compute the velocity of Freestream in the geometry axis system.
"""
velocity(fs::Freestream, ::Geometry) = velocity(fs)
# velocity(fs :: Freestream, ::Body) = flip_xz(velocity(fs, Geometry()))
# velocity(fs :: Freestream, ::Stability) = 
# velocity(fs :: Freestream, ::Wind) = 

include("reference_frames.jl")

geometry_to_wind_axes(xyz, fs::Freestream) = geometry_to_wind_axes(xyz, fs.alpha, fs.beta)
geometry_to_wind_axes(vor::AbstractVortex, fs::Freestream) =
    geometry_to_wind_axes(vor, fs.alpha, fs.beta)

function geometry_to_wind_axes(vortex::AbstractVortex, α, β)
    T = promote_type(eltype(α), eltype(β))
    return transform(vortex, LinearMap(RotZY{T}(β, α)))
end

function wind_to_geometry_axes(vor::AbstractVortex, α, β)
    T = promote_type(eltype(α), eltype(β))
    return transform(vor, LinearMap(RotYZ{T}(-α, -β)))
end

# Prandl-Glauert transformation
include("prandtl_glauert.jl")

## Matrix and residual setups
#==========================================================================================#

include("residuals.jl")

"""
    solve_linear(horseshoes, normals, U, Ω) 

Evaluate and return the vortex strengths ``Γ``s given an array of `Horseshoes`, their associated normal vectors, the velocity vector ``U``, and the quasi-steady rotation vector ``Ω``.
"""
function solve_linear(horseshoes, U, Ω)
    AIC  = influence_matrix(horseshoes)
    boco = boundary_condition(horseshoes, U, Ω)
    Γs   = AIC \ boco

    return Γs, AIC, boco
end

## Force evaluations
#==========================================================================================#

# Nearfield forces
include("nearfield.jl")

# Farfield forces
include("farfield.jl")

# System setups
#==========================================================================================#

# System
include("system.jl")

# Derivatives
include("stability.jl")

## Post-processing
#==========================================================================================#

# Pretty-printing
include("printing.jl")

# Streamlines
include("streamlines.jl")

end