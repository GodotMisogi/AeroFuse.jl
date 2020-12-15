module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs

## Math tools
#==========================================================================================#

include("../Tools/MathTools.jl")

using .MathTools: accumap, structtolist, three_quarter_point, quarter_point

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

## Vortex lattice
#==========================================================================================#

export solve_horseshoes

"""
    velocity(r, horseshoe, Γ, V_hat)

Computes the induced velocities at a point ``r`` of a given Horseshoe with constant strength ``Γ`` and trailing legs pointing in a given direction ``\\hat V``.
"""
velocity(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, Γ :: Real, V_hat :: SVector{3, <: Real}) = horseshoe_velocity(r, horseshoe.bound_leg, Γ, direction = V_hat)

"""
    mirror_velocity(r, horseshoe, Γ, V_hat)

Computes the induced velocity, using the method of images for a symmetric case in the x-z plane, at a point ``r`` of a given Horseshoe with constant strength ``Γ`` and trailing legs pointing in a given direction ``\\hat V``.
"""
function mirror_velocity(collocation_point :: SVector{3, <: Real}, horseshoe :: Horseshoe, Γ :: Real, V_hat :: SVector{3, <: Real})
    mirror_point = reflect_xz(collocation_point)
    mir_vel = (reflect_xz ∘ velocity)(mirror_point, horseshoe, Γ, V_hat)
end

"""
    influence_coefficient(collocation_point, horseshoe, panel_normal, V_hat, symmetry)

Computes the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V``` at a collocation point projected to a normal vector, with an option for symmetry in the x-z plane.
"""
function influence_coefficient(collocation_point :: SVector{3, <: Real}, horseshoe :: Horseshoe, panel_normal :: SVector{3, <: Real}, V_hat :: SVector{3, <: Real}, symmetry = false) 
    if symmetry
        col_vel = velocity(collocation_point, horseshoe, 1., V_hat)
        mir_vel = mirror_velocity(collocation_point, horseshoe, 1., V_hat)
        return dot(col_vel .+ mir_vel, panel_normal)
    else 
        return dot(velocity(collocation_point, horseshoe, 1., V_hat), panel_normal)
    end
end

# Matrix setup and solution
#==========================================================================================#

"""
    influence_matrix(colpoints, normals, horseshoes, V_hat, symmetry)

Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, a unit vector representing the freestream, with an option for symmetry in the x-z plane.
"""
influence_matrix(colpoints, normals, horseshoes :: AbstractVector{Horseshoe}, V_hat :: SVector{3, <: Real}, symmetry = false) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(colpoint_i, horsie_j, normal_i, V_hat, symmetry) for (colpoint_i, normal_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]

"""
    boundary_condition(velocities, normals)

Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = @timeit "Dotting" dot.(velocities, normals)

"""
    solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry) 

Solves the AIC matrix with the boundary condition given Panel3Ds and a Freestream, with the option to use the symmetry of the problem in the x-z plane.
"""
function solve_horseshoes(horseshoe_panels :: AbstractVector{Panel3D}, camber_panels :: AbstractVector{Panel3D}, freestream :: Freestream, symmetry = false) 
    U = aircraft_velocity(freestream)
    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = [ U .+ freestream.Ω × rc_i for rc_i ∈ colpoints ]

    @timeit "AIC" AIC = influence_matrix(colpoints[:], normals[:], horseshoes[:], -normalize(U), symmetry)
    @timeit "RHS" boco = boundary_condition(total_vel[:], normals[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = reshape(Γs, size(horseshoes)...), horseshoes
end


## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics

"""
Placeholder. Unsure whether to change this to a generic moment computation function.
"""
moments(horseshoes :: AbstractVector{Horseshoe}, forces, r_ref) = [ (bound_leg_center(vortex_ring) .- r_ref) × force for (force, vortex_ring) ∈ zip(forces, horseshoes) ]

## Post-processing
#==========================================================================================#

include("streamlines.jl")

export streamlines

end