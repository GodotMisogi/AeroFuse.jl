module AeroMDAO

#----------------------IMPORTS--------------------------------#
using StaticArrays
using Rotations

using TimerOutputs

## Math tools
include("General/math_tools.jl")

## Non-dimensionalization
#==========================================================================================#

include("General/nondimensional.jl")

export dynamic_pressure, force_coefficient, moment_coefficient, rate_coefficient, pressure_coefficient, aerodynamic_coefficients, print_dynamics

#-------------------------AIRCRAFT GUFF---------------------#

abstract type Aircraft end

## Foil geometry
#==========================================================================================#

include("Wings/Foils.jl")

export Foil, read_foil, kulfan_CST, naca4

## Wing geometry
#==========================================================================================#

include("Wings/wing_geometry.jl")

export HalfWing, Wing, mean_aerodynamic_chord, span, aspect_ratio, projected_area, info, print_info, lead_wing, wing_bounds, paneller, mesh_horseshoes, mesh_wing, mesh_cambers, make_panels, vlmesh_wing

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