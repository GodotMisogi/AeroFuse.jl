module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using ComponentArrays

## Package imports
#==========================================================================================#

# Math tools
import ..MathTools: weighted_vector, reshape_array

# Panel geometry
import ..PanelGeometry: Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform, p1, p2, p3, p4, average_chord, average_width

# Non-dimensionalization
import ..NonDimensional: dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient

# Tools regarding solutions to Laplace's equation
import ..Laplace: AbstractFreestream, Freestream, aircraft_velocity,cartesian_to_freestream, freestream_to_cartesian

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")
include("finite_core.jl")

## Reference frames
#==========================================================================================#

include("reference_frames.jl")

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

quasi_steady_freestream(horseshoes, U, Ω) = map(hs -> U + Ω × horseshoe_point(hs), horseshoes)

"""
    solve_system(horseshoes, normals, U, Ω) 

Evaluate and return the vortex strengths ``\\Gamma``s given `Horseshoes`, their associated normal vectors (not necessarily the same as the panels' normals), the speed ``U`` and rotation vector ``\\Omega``.
"""
function solve_system(horseshoes, U, Ω, finite_core = false)
    AIC  = influence_matrix(horseshoes, -normalize(U), finite_core)
    boco = boundary_condition(quasi_steady_freestream(horseshoes, U, Ω), horseshoe_normal.(horseshoes))
    Γs   = AIC \ boco 
end

## Force evaluations
#==========================================================================================#

# Nearfield forces
include("nearfield.jl")

# Farfield forces
include("farfield.jl")

## Post-processing
#==========================================================================================#

# Streamlines
include("streamlines.jl")

# System setups
#==========================================================================================#

# include("mutating_system.jl")
include("system.jl")

end