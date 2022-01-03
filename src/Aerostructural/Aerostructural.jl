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
import ..VortexLattice: velocity, induced_velocity, induced_trailing_velocity, aerodynamic_residuals, aerodynamic_residuals!, Horseshoe, horseshoe_normal, horseshoe_point, bound_leg_center, influence_coefficient, influence_matrix, boundary_condition, geometry_to_wind_axes, bound_leg_vector, kutta_joukowsky, surface_forces


struct AerostructWing{S,T}
    aerodynamics :: WingMesh{S}
    structures   :: Beam{T}
end

struct AerostructSystem{T}
    vortices :: Array{Horseshoe{T}}

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