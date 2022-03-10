module ViscFoil

import ..MathTools: midpair_map, forward_difference, forward_sum, weighted_vector

import ..Laplace: Uniform2D, velocity

import ..AeroMDAO: reynolds_number

import ..PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, distance, panel_length, panel_angle, tangent_vector, panel_points, panel_location, collocation_point

import ..DoubletSource: doublet_matrix, source_matrix, source_strengths, boundary_vector, solve_linear, lift_coefficient

import ..LinearVortexSource: vortex_influence_matrix, source_influence_matrix, linear_vortex_neumann_matrix, linear_source_neumann_matrix, neumann_boundary_condition

using NLsolve
using LineSearches: BackTracking
using ForwardDiff
using LinearAlgebra

include("inviscid.jl")

include("thermodynamics.jl")

include("closure_relations.jl")

include("blsys.jl")

include("viscous.jl")


include("system.jl")

end