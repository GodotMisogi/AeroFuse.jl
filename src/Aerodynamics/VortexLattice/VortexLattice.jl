module VortexLattice

using LinearAlgebra
using StaticArrays
using Rotations

## Panel geometry and math tools
#==========================================================================================#

import ..AeroMDAO: Panel3D, panel_area, panel_coords, midpoint, panel_normal, transform, point1, point2, point3, point4, dynamic_pressure, aerodynamic_coefficients, force_coefficient, moment_coefficient, rate_coefficient, fwddiff, weighted_vector

## Horseshoe setup
#==========================================================================================#

include("horseshoes.jl")

export Horseshoe, horseshoe_line, horseshoe_point

## Reference frames
#==========================================================================================#

include("reference_frames.jl")

export body_to_stability_axes, body_to_wind_axes, reflect_xz

## Influence matrix and solution of system
#==========================================================================================#

include("influences.jl")

export solve_horseshoes

"""
	make_horseshoes(horseshoe_panels)

Get the bound legs and collocation points of the horseshoes defined by `Panel3D`s.
"""
make_horseshoes(horseshoe_panels) =	@. horseshoe_line(horseshoe_panels), horseshoe_point(horseshoe_panels)

"""
	solve_horseshoes(horseshoe_panels, normals, U, Ω) 

Evaluate and return the vortex strengths ``\\Gamma``s and associated horsehoes given `Panel3D`s, their associated normal vectors (not necessarily the same as the panels' normals), the speed ``U`` and rotation vector ``\\Omega``.
"""
function solve_horseshoes(horseshoe_panels, normals, U, Ω) 
	horseshoes, colpoints = make_horseshoes(horseshoe_panels)
	total_vel = [ U + Ω × colpoint for colpoint in colpoints ]
	
	# Solve system
	Γs = solve_system(colpoints, horseshoes, normals, total_vel, U) 

	Γs, horseshoes
end

function solve_system(colpoints, horseshoes, normals, total_vel, U)
	AIC  = influence_matrix(colpoints, normals, horseshoes, -normalize(U))
	boco = boundary_condition(total_vel, normals)
	Γs 	 = AIC \ boco 
end

## Force evaluations
#==========================================================================================#

include("nearfield.jl")

export nearfield_dynamics, nearfield_drag

include("farfield.jl")

export farfield_dynamics

## Post-processing
#==========================================================================================#

export case_dynamics, evaluate_coefficients

"""
	case_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

Compute the nearfield and farfield forces over a given component with associated vortex strengths `Γ_comp` and horseshoes `hs_comp`, using the vortex strengths ``\\Gamma``s and `horseshoes` of the entire configuration, the speed ``U``, angles of attack ``\\alpha`` and sideslip ``\\beta``, the rotation vector ``\\Omega``, the freestream density ``\\rho`` and the reference location for moments ``r_\\text{ref}``.
"""
function case_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)
	# Compute near-field dynamics
	geom_forces, geom_moments = nearfield_dynamics(Γ_comp, hs_comp, Γs, horseshoes, U, Ω, ρ, r_ref)

	# Compute farfield dynamics
	trefftz_force = trefftz_forces(Γ_comp, hs_comp, norm(U), α, β, ρ)

	geom_forces, geom_moments, trefftz_force
end

# In case of only one component
case_dynamics(Γs, horseshoes, U, α, β, Ω, ρ, r_ref) = case_dynamics(Γs, horseshoes, Γs, horseshoes, U, α, β, Ω, ρ, r_ref)

function evaluate_coefficients(forces, moments, trefftz_force, U, α, β, Ω, ρ, S, c, b)
	V = norm(U)
	q = dynamic_pressure(ρ, V)

	# Computing summed coefficients
	force, moment = sum(forces), sum(moments)

	# Transform near-field dynamics to wind axes
	trans_force	 = body_to_wind_axes(force, α, β)
	trans_force  = [ nearfield_drag(force, U); trans_force[2:end] ]
	trans_moment = body_to_wind_axes(stability_flip(moment), α, β)

	# Compute coefficients
	nearfield_coeffs = aerodynamic_coefficients(trans_force, trans_moment, Ω, V, S, b, c, ρ)
	farfield_coeffs	 = force_coefficient(trefftz_force, q, S)

	# Non-dimensional panel coefficients
	CFs = force_coefficient.(forces, q, S)
	CMs = moment_coefficient.(moments, q, S, b, c)

	nearfield_coeffs, farfield_coeffs, CFs, CMs
end

# Streamlines
include("streamlines.jl")

export streamlines

end