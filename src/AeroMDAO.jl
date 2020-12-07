module AeroMDAO

#----------------------IMPORTS--------------------------------#
using StaticArrays
using Rotations

using TimerOutputs

## Math tools
#==========================================================================================#

include("General/math_tools.jl")
using .MathTools: tupvector, fwdsum, fwddiv, cosine_dist, weighted_vector, vectarray, slope, splitat, adj3, cosine_interp, columns

export tupvector


## Non-dimensionalization
#==========================================================================================#

include("General/nondimensional.jl")

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics

## Wing geometry
#==========================================================================================#

include("Wings/AircraftGeometry.jl")
# using .AircraftGeometry

## Vortex lattice
#==========================================================================================#

include("VortexLattice/VortexLattice.jl")
using .VortexLattice

export Panel3D, Freestream, velocity, streamlines, solve_horseshoes, transform

## Doublet-source
#==========================================================================================#

# include("DoubletSource/DoubletSource.jl")
# using .DoubletSource

# export Panel2D, Uniform2D, velocity

## Aerodynamic analyses
#==========================================================================================#

include("cases.jl")

export solve_case

## Post-processing
#==========================================================================================#

include("General/plot_tools.jl")

export plot_panels, plot_streamlines

end