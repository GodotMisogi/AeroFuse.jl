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

export body_to_stability_axes, body_to_wind_axes, reflect_xz

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

export solve_horseshoes

"""
Placeholder.
"""
function make_horseshoes(horseshoe_panels :: AbstractVector{Panel3D})
    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)

    horseshoes, colpoints
end

function solve_system(colpoints, horseshoes, normals, total_vel, U, symmetry)
    # Allocated versions
    # @timeit "AIC" AIC = influence_matrix(colpoints, normals, horseshoes, -normalize(U), symmetry)
    # @timeit "RHS" boco = boundary_condition(total_vel, normals)
    
    # Pre-allocated versions
    AIC = zeros(length(colpoints), length(colpoints))
    boco = zeros(length(colpoints))

    @timeit "AIC (Preallocated)" influence_matrix!(AIC, colpoints, normals, horseshoes, -normalize(U))
    @timeit "RHS (Preallocated)" boundary_condition!(boco, total_vel, normals)

    @timeit "Solve AIC" AIC \ boco
end

"""
    solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry) 

Solves the AIC matrix with the boundary condition given Panel3Ds and a Freestream, with the option to use the symmetry of the problem in the ``x``-``z`` plane.
"""
function solve_horseshoes(horseshoe_panels :: AbstractVector{Panel3D}, camber_panels :: AbstractVector{Panel3D}, freestream :: Freestream, symmetry = false) 
    @timeit "Freestream Velocity" U = aircraft_velocity(freestream)

    horseshoes, colpoints = make_horseshoes(horseshoe_panels)
    
    @timeit "Normals" normals = panel_normal.(camber_panels)
    
    @timeit "Total Velocity" total_vel = [ U + freestream.Ω × colpoint for colpoint in colpoints ]

    @timeit "Solving..." Γs = solve_system(colpoints, horseshoes, normals, total_vel, U, symmetry)

    @timeit "Reshape" Γs, horseshoes
end


## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics

## Post-processing
#==========================================================================================#

export case_dynamics

function case_dynamics(Γs :: AbstractArray{<: Real}, horseshoes :: AbstractArray{Horseshoe}, freestream :: Freestream, r_ref :: SVector{3, <: Real}, ρ :: Real, symmetry :: Bool = false)
    # Compute near-field dynamics
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs[:], horseshoes[:], freestream, r_ref, ρ, symmetry)
    force, moment = sum(geom_forces), sum(geom_moments)
    drag = nearfield_drag(force, freestream)

    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream)

    # Compute farfield dynamics
    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ, symmetry)

    geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment
end

include("streamlines.jl")

export streamlines

end