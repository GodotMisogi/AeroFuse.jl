module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations, CoordinateTransformations
using ComponentArrays
using SplitApplyCombine
using TimerOutputs
using LabelledArrays
using Setfield

## Package imports
#==========================================================================================#

# Math tools
import ..MathTools: weighted_vector, structtolist

# Panel geometry
import ..PanelGeometry: Panel3D, panel_area, panel_coordinates, midpoint, panel_normal, transform, p1, p2, p3, p4, average_chord, average_width

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Some tools
import ..Laplace: cartesian_to_freestream, freestream_to_cartesian

import ..AeroMDAO: velocity, solve_system, solve_linear, surface_velocities, surface_coefficients

## Vortex types and methods
#==========================================================================================#

include("horseshoes.jl")

include("vortex_rings.jl")

## Reference frames
#==========================================================================================#

include("freestream.jl")

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
    AIC  = influence_matrix(horseshoes, -normalize(U))
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

include("system.jl")

function solve_system(components, fs :: Freestream, refs :: References)

    # Mach number bound checks
    M = mach_number(refs)
    @assert M < 1.  "Only compressible subsonic flow conditions (M < 1) are valid!"
    if M > 0.7 @warn "Results in transonic flow conditions (0.7 < M < 1) are most likely incorrect!" end

    # # Prandtl-Glauert transformation
    # β          = √(1 - M^2)
    # components = @. prandtl_glauert_scale_coordinates(geometry_to_wind_axes(components, -fs.alpha, -fs.beta), β)
    # U          = wind_to_geometry_axes(body_frame_velocity(fs), -fs.alpha, -fs.beta)
    # Ω          = wind_to_geometry_axes(fs.omega, -fs.alpha, -fs.beta)

    U, Ω = body_frame_velocity(fs), fs.omega

    # Solve system
    Γs, AIC, boco = solve_linear(components, U, Ω)

    VortexLatticeSystem(components, refs.speed * Γs, AIC, boco, fs, refs)
end

## Post-processing
#==========================================================================================#

# Streamlines
include("streamlines.jl")


end