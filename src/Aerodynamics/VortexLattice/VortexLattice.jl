module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations, CoordinateTransformations
using ComponentArrays
using SplitApplyCombine
using TimerOutputs
using LabelledArrays
using Setfield
using PrettyTables
using ForwardDiff: jacobian!
using DiffResults: JacobianResult, jacobian, value

## Package imports
#==========================================================================================#

# Math tools
import ..MathTools: weighted_vector, structtolist, Point3D

# Panel geometry
import ..PanelGeometry: Panel3D, panel_area, panel_coordinates, midpoint, normal_vector, transform, p1, p2, p3, p4, average_chord, average_width

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Some tools
import ..Laplace: AbstractFreestream, Freestream

import ..AeroFuse: velocity, solve_system, solve_linear, solve_nonlinear, solve_nonlinear!, surface_velocities, surface_coefficients

## Vortex types and methods
#==========================================================================================#

include("horseshoes.jl")

include("vortex_rings.jl")

## Reference frames
#==========================================================================================#

## Axis transformations
abstract type AbstractAxisSystem end

struct Geometry <: AbstractAxisSystem end
struct Body      <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind      <: AbstractAxisSystem end

Base.show(io :: IO, :: Geometry)  = print(io, "Geometry")
Base.show(io :: IO, :: Body)      = print(io, "Body")
Base.show(io :: IO, :: Stability) = print(io, "Stability")
Base.show(io :: IO, :: Wind)      = print(io, "Wind")

include("freestream.jl")

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

# CONVERT THIS TO VLM CONSTRUCTOR
function solve_system(components, fs :: Freestream, refs :: References)

    # Mach number bound checks
    M = mach_number(refs)
    @assert M < 1.  "Only compressible subsonic flow conditions (M < 1) are valid!"
    if M > 0.7 @warn "Results in transonic flow conditions (0.7 < M < 1) are most likely incorrect!" end

    # Wind axis + Prandtl-Glauert transformation
    β_pg = √(1 - M^2)
    comp = @. prandtl_glauert_scale_coordinates(geometry_to_wind_axes(components, fs), β_pg)

    # Quasi-steady freestream velocity
    U = geometry_to_wind_axes(velocity(fs, Body()), fs)
    Ω = geometry_to_wind_axes(fs.omega, fs) / refs.speed

    # Solve system
    Γs, AIC, boco = solve_linear(comp, U, Ω)

    return VortexLatticeSystem(components, refs.speed * Γs / β_pg^2, AIC, boco, fs, refs)
end

## Post-processing
#==========================================================================================#

# Pretty-printing
include("printing.jl")

# Streamlines
include("streamlines.jl")


end