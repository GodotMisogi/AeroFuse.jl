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

import ..AeroFuse: velocity, solve_system, solve_linear, solve_nonlinear, solve_nonlinear!, surface_velocities, surface_coefficients

## Vortex types and methods
#==========================================================================================#

include("vortices.jl")

## Reference frames
#==========================================================================================#

## Axis transformations with respect to freestream
abstract type AbstractAxisSystem end

struct Geometry <: AbstractAxisSystem end
struct Body      <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind      <: AbstractAxisSystem end

"""
    velocity(freestream :: Freestream, ::Body)

Compute the velocity of Freestream in the body reference frame.
"""
velocity(fs :: Freestream, ::Body) = -velocity(fs)

include("reference_frames.jl")

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

    Γs, AIC, boco
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