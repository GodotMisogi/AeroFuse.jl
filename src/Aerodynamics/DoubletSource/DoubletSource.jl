module DoubletSource

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics

using ..AeroMDAO: Panel, Panel2D, WakePanel2D, Point2D, collocation_point, point1, point2, transform_panel, affine_2D, panel_length, panel_angle, panel_tangent, panel_normal, panel_dist, rotation, inverse_rotation, midpair_map, pressure_coefficient, wake_panel, wake_panels, panel_points, panel_vector

## Doublet-source Dirichlet boundary condition
#===========================================================================#

source_potential(str, x, z, x1, x2) = str / 4π * ((x - x1) * log((x - x1)^2 + z^2) - (x - x2) * log((x - x2)^2 + z^2) + 2z * (atan(z, x - x2) - atan(z, x - x1)))

source_velocity(str, x, z, x1, x2) = SVector(str / (4π) * log(((x - x1)^2 + z^2) / ((x - x2)^2 + z^2)), doublet_potential(str, x, z, x1, x2))

doublet_potential(str, x, z, x1, x2) = str / (2π) * (atan(z, x - x1) - atan(z, x - x2))

doublet_velocity(str, x, z, x1, x2) = SVector(str / (2π) * - (z / ((x - x1)^2 + z^2) - z / ((x - x2)^2 + z^2) ), str / (2π) * ( (x - x1) / ((x - x1)^2 + z^2) - (x - x2) / ((x - x2)^2 + z^2)))


## Matrix helpers
#===========================================================================#

function doublet_influence(panel_j :: Panel, panel_i :: Panel)
	xp, yp = transform_panel(panel_j, panel_i)
	doublet_potential(1., xp, yp, 0., panel_length(panel_j))
end

function source_influence(panel_j :: Panel, panel_i :: Panel)
	xp, yp = transform_panel(panel_j, panel_i)
	source_potential(1., xp, yp, 0., panel_length(panel_j))
end

boundary_condition(panel_j :: Panel, panel_i :: Panel, u) = -source_influence(panel_j, panel_i) * dot(u, panel_normal(panel_j))

## Aerodynamic coefficients
#===========================================================================#

panel_velocity(dφ, dr, u, α) = dφ / dr + dot(u, α)

lift_coefficient(cp, dist_colpoints, panel_angle) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength, speed) = 2. * wake_strength / speed

"""
	aerodynamic_coefficients(vels, Δrs, panel_angles, speed, α)

Compute the lift, moment, and pressure coefficients given associated arrays of edge speeds, adjacent collocation point distances, panel angles, the freestream speed, and angle of attack ``\\alpha``.
"""
function eval_coefficients(vels, Δrs, xjs, panel_angles, speed, α)
	cps   = @. pressure_coefficient(speed, vels)
	cls   = @. lift_coefficient(cps, Δrs, panel_angles)
	cms   = @. -cls * xjs * cos(α)

	cls, cms, cps
end

## Matrix assembly
#===========================================================================#

include("matrix_func.jl")

export solve_problem

function solve_problem(panels :: Vector{<: Panel2D}, u, α, sources :: Bool, wake_length)
	speed 			= norm(u)
	xs	 			= getindex.(panel_points(panels)[2:end-1], 1)

	# Blunt trailing edge tests
	te_panel 		= Panel2D((point2 ∘ last)(panels), (point1 ∘ first)(panels))
	r_te 			= panel_vector(te_panel)
	φ_TE  			= dot(u, r_te)

	# Solve for doublet strengths
	φs 				= solve_strengths(panels, u, α, r_te, sources; bound = wake_length)
	
	# Evaluate inviscid edge velocitiess
	u_es, Δrs 		= tangential_velocities(panels, φs, u, sources)
	
	# Compute coefficients 
	cls, cms, cps 	= eval_coefficients(u_es, Δrs, xs, panel_angle.(panels[2:end]), speed, α)
	
	# Evaluate lift coefficient from wake doublet strength
	cl_wake 		= lift_coefficient(φs[end] - φs[1] + φ_TE, speed)

	cls, cms, cps, cl_wake
end

function solve_problem(panels :: Vector{<: Panel2D}, u, α, num_wake :: Integer, wake_length)
	speed 			= norm(u)
	xs	 			= getindex.(panel_points(panels)[2:end-1], 1)
	wakes 			= wake_panels(panels, wake_length, num_wake)
	φs				= solve_strengths(panels, u, α, wakes; bound = wake_length)
	u_es, Δrs 		= tangential_velocities(panels, φs[1:end-1], u, false)
	cls, cms, cps 	= evaluate_coefficients(u_es, Δrs, xs, panel_angle.(panels[2:end]), speed, α)
	cl_wake 		= lift_coefficient(last(φs), speed)

	cls, cms, cps, cl_wake
end


end