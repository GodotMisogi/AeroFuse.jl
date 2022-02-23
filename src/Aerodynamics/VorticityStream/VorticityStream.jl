module VorticityStream

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using TimerOutputs

import ..PanelGeometry: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, collocation_point, p1, p2, transform_panel, transform_panel_points, affine_2D, panel_length, panel_angle, panel_vector, tangent_vector, normal_vector, distance, panel_points, reverse_panel, panel_scalar, trailing_edge_panel

import ..Laplace: Uniform2D

import ..NonDimensional: pressure_coefficient

import ..MathTools: rotation, inverse_rotation, midpair_map

import ..AeroMDAO: solve_linear, solve_linear!

include("singularities.jl")
include("matrix.jl")

struct VorticityStreamSystem{N,T}
    influence_matrix :: SMatrix{N,N,T}
    boundary_vector  :: SMatrix{N,2,T}
    singularities    :: SMatrix{N,2,T}
end

# export solve_problem

# function solve_problem(panels :: Vector{<: Panel2D}, u, sources :: Bool, wake_length)
#   φs        = solve_linear(panels, u, sources; bound = wake_length)
#   cps, cls  = aerodynamic_coefficients(panels, φs, u, sources)
#   cl_wake   = lift_coefficient(last(φs), norm(u))

#   cps, cls, cl_wake
# end

# function solve_problem(panels :: Vector{<: Panel2D}, u, num_wake :: Integer, wake_length)
#   wakes     = wake_panels(panels, wake_length, num_wake)
#   φs    = solve_linear(panels, u, wakes; bound = wake_length)
#   cps, cls  = aerodynamic_coefficients(panels, φs, u, true)
#   cl_wake   = lift_coefficient(last(φs), norm(u))

#   cps, cls, cl_wake
# end



# function ψp_influence(panel_i :: Panel2D, point)
#   xp, yp = transform_panel(panel_i, point)
#   vortex_stream_1(1., xp, yp, 0., panel_length(panel_i)) + vortex_stream_2(1., xp, yp, 0., panel_length(panel_i))
# end

# function ψm_influence(panel_i :: Panel2D, point)
#   xp, yp = transform_panel(panel_i, point)
#   vortex_stream_1(1., xp, yp, 0., panel_length(panel_i)) - vortex_stream_2(1., xp, yp, 0., panel_length(panel_i))
# end

# function te_influence(panel_i :: Panel2D, point)
#   xp, yp = transform_panel(panel_i, point)
#   source_stream(1., xp, yp, 0., panel_length(panel_i)) + vortex_stream_1(1., xp, yp, 0., panel_length(panel_i))
# end

# function source_influence(panel_i :: Panel2D, point)
#   xp, yp = transform_panel(panel_i, point)
#   source_stream(1., xp, yp, 0., panel_length(panel_i))
# end

# function source_influence_te(panel_i :: Panel2D, point)
#   xp, yp = transform_panel(panel_i, point)
#   source_stream_te(1., xp, yp, 0., panel_length(panel_i))
# end

end