module AeroMDAO

#----------------------IMPORTS--------------------------------#
using StaticArrays
using Rotations
using LinearAlgebra

using TimerOutputs

## Math tools
#==========================================================================================#

include("Tools/MathTools.jl")
using .MathTools: tupvector, fwdsum, fwddiv, cosine_dist, weighted_vector, vectarray, slope, splitat, adj3, cosine_interp, columns


## Non-dimensionalization
#==========================================================================================#

include("Tools/NonDimensional.jl")
using .NonDimensional

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics, reynolds_number

## Wing geometry
#==========================================================================================#

include("Geometry/AircraftGeometry.jl")
# using .AircraftGeometry

## Vortex lattice
#==========================================================================================#

include("VortexLattice/VortexLattice.jl")
using .VortexLattice

export Panel3D, Horseshoe, Freestream, velocity, streamlines, solve_horseshoes, transform, panel_coords

## Doublet-source panel method
#==========================================================================================#

include("DoubletSource/DoubletSource.jl")
using .DoubletSource

export Panel2D, Uniform2D, lift_coefficient

## Aerodynamic analyses
#==========================================================================================#

include("cases.jl")

export solve_case

## Post-processing
#==========================================================================================#

include("Tools/plot_tools.jl")

export plot_panels, plot_streamlines, trace_surface, trace_panels, trace_coords, trace_streamlines, panel_splits

end