module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations
using TimerOutputs

## Math tools
#==========================================================================================#

include("../General/MathTools.jl")

using .MathTools: accumap, structtolist, three_quarter_point, quarter_point

## Freestream
#==========================================================================================#

include("../General/laplace.jl")

export Laplace, Freestream, velocity, aircraft_velocity

## Reference frames
#==========================================================================================#

include("../General/reference_frames.jl")

export body_to_stability_axes, body_to_wind_axes

## Panel geometry
#==========================================================================================#

include("../General/PanelGeometry.jl")

import .PanelGeometry: Panel, Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform

export Panel, Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform

## Vortex lattice
#==========================================================================================#

export solve_horseshoes

"""
Helper function to compute the bound leg of Panel3D for horseshoes/vortex rings, which is the quarter point on each side in the x-z plane.
"""
bound_leg(p1, p2, p3, p4) = [ quarter_point(p1, p2), quarter_point(p4, p3) ]

"""
Helper function to compute the collocation point of Panel3D for horseshoes/vortex rings, which is the 3-quarter point on each side in the x-z plane.
"""
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) .+ three_quarter_point(p4, p3) ) ./ 2

"""
A composite type consisting of two vectors to define a line.
"""
struct Line
    r1 :: SVector{3, Float64}
    r2 :: SVector{3, Float64}
end

vector(line :: Line) = line.r2 - line.r1
center(line :: Line) = (line.r1 .+ line.r2) / 2

transform(line :: Line, rotation, translation) = Line(rotation * line.r1 + translation, rotation * line.r2 + translation)

"""
Helper function to compute the velocity induced by a bound vortex leg.
"""
bound_leg_velocity(a, b, Γ) = Γ/4π * (1/norm(a) + 1/norm(b)) * cross(a, b) / (norm(a) * norm(b) + dot(a, b))

"""
Helper function to compute the velocity induced by trailing vortex legs.
"""
trailing_legs_velocities(a, b, Γ, û) = Γ/4π * (cross(a, û) / (norm(a) - dot(a, û)) / norm(a) - cross(b, û) / (norm(b) - dot(b, û)) / norm(b))

"""
Helper function to check if any point is on the line.
"""
on_line(a, b, ε) = any(<(ε), norm.([a, b, cross(a, b) ]))

"""
Computes the velocity induced at a point `r` by a vortex Line with constant strength Γ. Checks if `r` lies on the line and sets it to `(0, 0, 0)` if so, as the velocity is singular there.
"""
function horseshoe_velocity(r, line :: Line, Γ, ε = 1e-8; direction = SVector(1, 0, 0))
    a, b = r .- line.r1, r .- line.r2

    # Compute bound leg velocity
    @timeit "Bound Leg" bound_velocity = on_line(a, b, ε) ? zeros(3) : bound_leg_velocity(a, b, Γ)
    
    # Compute velocities of trailing legs
    @timeit "Trailing Leg" trailing_velocity = trailing_legs_velocities(a, b, Γ, direction)

    # Sums bound and trailing legs' velocities
    @timeit "Sum Legs" bound_velocity .+ trailing_velocity
end

"""
Placeholder. Vortex rings and horseshoes basically have the same methods, and are arrays of vortex lines.
"""
abstract type AbstractVortexArray end

"""
A horseshoe type consisting of a bound leg of type Line represening a vortex line.
"""
struct Horseshoe <: AbstractVortexArray
    bound_leg :: Line
end

"""
Computes the bound leg for a Panel3D.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Computes the collocation point of a Panel3D.
"""
collocation_point(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Returns a Horseshoe on a Panel3D.
"""
horseshoe_lines(panel :: Panel3D) = (Horseshoe ∘ Line)(bound_leg(panel)...)

"""
Computes the midpoint of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_center(vortex :: AbstractVortexArray) = center(vortex.bound_leg)

"""
Computes the direction vector of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_vector(vortex :: AbstractVortexArray) = vector(vortex.bound_leg)


# Matrix setup
#==========================================================================================#

"""
Computes the velocity using the method of images for a symmetric case in the x-z plane.
"""
function mirror_velocity(collocation_point, horseshoe, Γ, Û)
    mirror_point = reflect_xz(collocation_point)
    mir_vel = (reflect_xz ∘ velocity)(mirror_point, horseshoe, Γ, Û)
end

"""
Computes the induced velocities at a point `r` of a Horseshoe with constant strength Γ and trailing legs pointing in a given direction.
"""
velocity(r, horseshoe :: Horseshoe, Γ, V̂) = horseshoe_velocity(r, horseshoe.bound_leg, Γ, direction = V̂)

"""
Computes the influence coefficient of the velocity of a vortex line at a collocation point projected to a normal vector.
"""
function influence_coefficient(collocation_point, horseshoe :: Horseshoe, panel_normal, V̂, symmetry = false) 
    if symmetry
        col_vel = velocity(collocation_point, horseshoe, 1., V̂)
        mir_vel = mirror_velocity(collocation_point, horseshoe, 1., V̂)
        return dot(col_vel .+ mir_vel, panel_normal)
    else 
        return dot(velocity(collocation_point, horseshoe, 1., V̂), panel_normal)
    end
end

"""
Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points and associated normal vectors.
"""
influence_matrix(colpoints, normals, horseshoes :: Array{Horseshoe}, V̂, symmetry = false) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" influence_coefficient(colpoint_i, horsie_j, normal_i, V̂, symmetry) for (colpoint_i, normal_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]

"""
Computes the projection of a velocity vector with respect to normal vectors of panels. Corresponds to construction of the boundary condition for the RHS of the AIC system.
"""
boundary_condition(velocities, normals) = @timeit "Dotting" dot.(velocities, normals)

"""
Solves the AIC matrix with the boundary condition given Panel3Ds and a freestream velocity.
"""
function solve_horseshoes(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, freestream :: Freestream, Ω, symmetry = false) 
    U = aircraft_velocity(freestream)
    @timeit "Horseshoes" horseshoes = horseshoe_lines.(horseshoe_panels)
    @timeit "Collocation Points" colpoints = collocation_point.(horseshoe_panels)
    @timeit "Normals" normals = panel_normal.(camber_panels)
    @timeit "Total Velocity" total_vel = [ U .+ cross(freestream.Ω, rc_i) for rc_i ∈ colpoints ]

    @timeit "AIC" AIC = influence_matrix(colpoints[:], normals[:], horseshoes[:], -normalize(U), symmetry)
    @timeit "RHS" boco = boundary_condition(normals[:], total_vel[:])
    @timeit "Solve AIC" Γs = AIC \ boco

    @timeit "Reshape" output = reshape(Γs, size(horseshoes)...), horseshoes

    output
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
moments(horseshoes :: Array{Horseshoe}, forces, r_ref) = [ cross((bound_leg_center(vortex_ring) .- r_ref), force) for (force, vortex_ring) ∈ zip(forces, horseshoes) ]

## Post-processing
#==========================================================================================#

include("streamlines.jl")

export streamlines

end