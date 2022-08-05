module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using ComponentArrays
using SplitApplyCombine
using TimerOutputs
using LabelledArrays

## Package imports
#==========================================================================================#

# Math tools
import ..MathTools: weighted_vector, structtolist

# Panel geometry
import ..PanelGeometry: Panel3D, panel_area, panel_coordinates, midpoint, panel_normal, p1, p2, p3, p4, average_chord, average_width

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

include("lines.jl")
include("reference_frames.jl")

## Matrix and residual setups
#==========================================================================================#

include("residuals.jl")

"""
    solve_linear(horseshoes, normals, U, Ω) 

Evaluate and return the vortex strengths ``Γ``s given `Horseshoes`, their associated normal vectors (not necessarily the same as the panels' normals), the speed ``U`` and rotation vector ``Ω``.
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
    Γs, AIC, boco = solve_linear(components, body_frame_velocity(fs), fs.omega)

    VortexLatticeSystem(components, refs.speed * Γs, AIC, boco, fs, refs)
end

## Post-processing
#==========================================================================================#

# Streamlines
include("streamlines.jl")


end