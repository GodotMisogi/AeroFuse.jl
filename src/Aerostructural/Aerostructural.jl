module Aerostructural

## Package imports
#==========================================================================================#

using StaticArrays
using LinearAlgebra
using ComponentArrays
using TimerOutputs
using SparseArrays
import SplitApplyCombine: combinedimsview, splitdimsview, filterview

# Structures
import ..Beams: tube_stiffness_matrix, AbstractBeam, Beam, Tube

# Conversions
import ..Laplace: freestream_to_cartesian

# Panelling
import ..PanelGeometry: make_panels, panel_normal

import ..AircraftGeometry: WingMesh

# VLM Aerodynamics
import ..VortexLattice: velocity, trailing_velocity, Horseshoe, horseshoe_normal, horseshoe_point, bound_leg_center, quasi_steady_freestream, influence_coefficient, influence_matrix, boundary_condition, surface_forces, body_to_wind_axes, surface_velocity, bound_leg_vector, kutta_joukowsky


struct AerostructWing{S,T}
    aerodynamics :: WingMesh{S}
    structures   :: Beam{T}
end

## Aerodynamic analysis
#==========================================================================================#

include("aerodynamics.jl")

## Structural analysis
#==========================================================================================#

include("structures.jl")

## Load-displacement transfer mechanisms
#==========================================================================================#

include("transfers.jl")

## Weights, engine, and fuel loads
#==========================================================================================#

# ???

## Coupled residual systems
#==========================================================================================#

include("residuals.jl")

end