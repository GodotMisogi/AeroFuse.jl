module DoubletSource

## Package imports
#==========================================================================================#

using LinearAlgebra
import Base.Iterators: product
using StaticArrays
import SplitApplyCombine: combinedimsview

import ..MathTools: rotation, inverse_rotation, midpair_map, Point3D

import ..Laplace: Uniform2D, magnitude, angle, velocity

import ..NonDimensional: pressure_coefficient

import ..PanelGeometry: AbstractPanel2D, Panel2D, WakePanel2D, collocation_point, p1, p2, transform_panel, affine_2D, panel_length, panel_angle, panel_tangent, panel_normal, distance, wake_panel, wake_panels, panel_points, panel_vector, AbstractPanel3D, Panel3D, panel_coordinates, midpoint, xs, ys, zs

import ..AeroMDAO: solve_system, surface_velocities, surface_coefficients

include("singularities3D.jl")

## Doublet-source Dirichlet boundary condition
#===========================================================================#

source_potential(str, x, z, x1, x2) = str / 4π * ((x - x1) * log((x - x1)^2 + z^2) - (x - x2) * log((x - x2)^2 + z^2) + 2z * (atan(z, x - x2) - atan(z, x - x1)))

source_velocity(str, x, z, x1, x2) = SVector(str / (4π) * log(((x - x1)^2 + z^2) / ((x - x2)^2 + z^2)), doublet_potential(str, x, z, x1, x2))

doublet_potential(str, x, z, x1, x2) = str / (2π) * (atan(z, x - x1) - atan(z, x - x2))

doublet_velocity(str, x, z, x1, x2) = SVector(str / (2π) * - (z / ((x - x1)^2 + z^2) - z / ((x - x2)^2 + z^2) ), str / (2π) * ( (x - x1) / ((x - x1)^2 + z^2) - (x - x2) / ((x - x2)^2 + z^2)))


## Matrix helpers
#===========================================================================#

function doublet_influence(panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D)
    xp, yp = transform_panel(panel_j, panel_i)
    ifelse(panel_i == panel_j, 0.5, doublet_potential(1., xp, yp, 0., panel_length(panel_j)))
end

function doublet_influence(panel_j :: AbstractPanel3D, panel_i :: AbstractPanel3D)
    xp, yp = transform_panel(panel_j, panel_i)
    ifelse(panel_i == panel_j, 0.5, doublet_potential(1., xp, yp, 0., panel_length(panel_j)))
end

function source_influence(panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D)
    xp, yp = transform_panel(panel_j, panel_i)
    source_potential(1., xp, yp, 0., panel_length(panel_j))
end

boundary_condition(panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D, u) = -source_influence(panel_j, panel_i) * dot(u, panel_normal(panel_j))

## Aerodynamic coefficients
#===========================================================================#

surface_velocity(dφ, dr, u, α) = dφ / dr + dot(u, α)

# """
#     aerodynamic_coefficients(vels, Δrs, panel_angles, speed, α)

# Compute the lift, moment, and pressure coefficients given associated arrays of edge speeds, adjacent collocation point distances, panel angles, the freestream speed, and angle of attack ``α``.
# """
# function evaluate_coefficients(vels, Δrs, xjs, panel_angles, speed, α)
#     cps   = @. pressure_coefficient(speed, vels)
#     cls   = @. lift_coefficient(cps, Δrs, panel_angles)
#     cms   = @. -cls * xjs * cos(α)

#     cls, cms, cps
# end

## Matrix assembly
#===========================================================================#

include("matrix_func.jl")

struct DoubletSourceSystem{T <: Real, M <: AbstractMatrix{T}, N <: AbstractVector{T}, O <: AbstractVector{<: AbstractPanel2D}, R <: WakePanel2D, P <: Uniform2D}
    influence_matrix   :: M
    boundary_condition :: N
    singularities      :: N
    surface_panels     :: O
    wake_panels        :: R
    freestream         :: P
end

function Base.show(io :: IO, sys :: DoubletSourceSystem)
    println(io, "DoubletSourceSystem —")
    println(io, length(sys.surface_panels), " ", eltype(sys.surface_panels), " Elements")
end

function solve_system(panels, uni :: Uniform2D, sources :: Bool, wake_length)
    # Freestream conditions
    u, α  = velocity(uni), uni.angle

    # Build wake
    wake_pan = wake_panel(panels, wake_length, α)

    # speed           = norm(u)
    # xs              = getindex.(panel_points(panels)[2:end-1], 1)

    # Blunt trailing edge tests
    # te_panel        = Panel2D((p2 ∘ last)(panels), (p1 ∘ first)(panels))
    # r_te            = panel_vector(te_panel)
    # φ_TE            = dot(u, r_te)

    
    # Solve for doublet strengths
    φs, AIC, boco   = solve_linear(panels, u, α, wakes; bound = wake_length)

    DoubletSourceSystem(AIC, boco, φs, panels, wakes, uni)

    # # Evaluate inviscid edge velocities
    # u_es, Δrs       = tangential_velocities(panels, φs, u, sources)

    # # Compute coefficients
    # cls, cms, cps   = evaluate_coefficients(u_es, Δrs, xs, panel_angle.(panels[2:end]), speed, α)

    # # Evaluate lift coefficient from wake doublet strength
    # cl_wake         = lift_coefficient(φs[end] - φs[1] + φ_TE, speed)

    # cls, cms, cps, cl_wake
end

function solve_system(panels, uni :: Uniform2D, num_wake :: Integer, wake_length)
    u, α  = velocity(uni), uni.angle

    wake_pan = wake_panel(panels, wake_length, α)
    # wakes = wake_panels(panels, wake_length, num_wake)

    # Solve for doublet strengths
    φs, AIC, boco   = solve_linear(panels, u, wake_pan) # ; bound = wake_length)

    DoubletSourceSystem(AIC, boco, φs, panels, wake_pan, uni)
end

function surface_velocities(prob :: DoubletSourceSystem)
    # Panel properties
    ps   = prob.surface_panels
    Δrs  = @views @. distance(ps[2:end], ps[1:end-1])
    αs   = @views panel_tangent.(ps[2:end])
    
    @views surface_velocities(prob.singularities[1:end-1], Δrs, αs, velocity(prob.freestream), false)
end

function surface_coefficients(prob :: DoubletSourceSystem)
    # Panel properties
    ps   = prob.surface_panels
    Δrs  = @views @. distance(ps[2:end], ps[1:end-1])
    xs   = @views combinedimsview(panel_points(ps)[2:end-1])[1,:]
    θs   = @views panel_angle.(ps[2:end])

    # Inviscid edge velocities
    u_es = @views surface_velocities(prob)

    # Aerodynamic coefficients
    cps  = @. pressure_coefficient(prob.freestream.magnitude, u_es)
    cls  = @. -cps * Δrs * cos(θs)
    cms  = @. -cls * xs * cos(prob.freestream.angle)

    cls, cms, cps
end

lift_coefficient(prob :: DoubletSourceSystem) = 2 * last(prob.singularities) / prob.freestream.magnitude

end