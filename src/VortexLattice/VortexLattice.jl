module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs

## Math tools
#==========================================================================================#

include("../Tools/MathTools.jl")

using .MathTools: accumap, fwddiff, structtolist, three_quarter_point, quarter_point

## Freestream
#==========================================================================================#

include("../Tools/Laplace.jl")
import .Laplace: Freestream, velocity, aircraft_velocity

export Freestream, velocity, aircraft_velocity

## Panel geometry
#==========================================================================================#

include("../Geometry/PanelGeometry.jl")
import .PanelGeometry: Panel, Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform

export Panel, Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")

export Horseshoe

## Reference frames
#==========================================================================================#

include("../Tools/ReferenceFrames.jl")
# using .ReferenceFrames: body_to_stability_axes, body_to_wind_axes

export body_to_stability_axes, body_to_wind_axes

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")
include("symmetry.jl")

export solve_horseshoes

"""
    solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry) 

Solves the AIC matrix with the boundary condition given Panel3Ds and a Freestream, with the option to use the symmetry of the problem in the ``x``-``z`` plane.
"""
function solve_horseshoes(horseshoe_panels :: AbstractVector{Panel3D}, camber_panels :: AbstractVector{Panel3D}, freestream :: Freestream, symmetry = false) 
    U = aircraft_velocity(freestream)
    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = Ref(U) .+ Ref(freestream.Ω) .× colpoints
    
    @timeit "AIC" AIC = symmetry ? sym_influence_matrix(colpoints[:], normals[:], horseshoes[:], -normalize(U)) : influence_matrix(colpoints[:], normals[:], horseshoes[:], -normalize(U))
    @timeit "RHS" boco = boundary_condition(total_vel[:], normals[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = Γs, horseshoes
end


## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics


## Post-processing
#==========================================================================================#

include("streamlines.jl")

export streamlines

end