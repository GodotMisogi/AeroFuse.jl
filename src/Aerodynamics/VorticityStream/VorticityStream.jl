module VorticityStream

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using TimerOutputs

using ..AeroMDAO: AbstractPanel, AbstractPanel2D, Panel2D, WakePanel2D, Point2D, Uniform2D, collocation_point, p1, p2, transform_panel, transform_panel_points, affine_2D, panel_length, panel_angle, panel_tangent, panel_normal, panel_dist, panel_points, stream, pressure_coefficient, rotation, inverse_rotation, midpair_map, panel_points, reverse_panel, panel_scalar, trailing_edge_panel

include("singularities.jl")
include("matrix.jl")

## Dynamics helpers
#===========================================================================#

lift_coefficient(cp, dist_colpoints, panel_angle) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength, speed) = 2. * wake_strength / speed

## Matrix assembly
#===========================================================================#


# export solve_problem

# function solve_problem(panels :: Vector{<: Panel2D}, u, sources :: Bool, wake_length)
#   φs        = solve_strengths(panels, u, sources; bound = wake_length)
#   cps, cls  = aerodynamic_coefficients(panels, φs, u, sources)
#   cl_wake   = lift_coefficient(last(φs), norm(u))

#   cps, cls, cl_wake
# end

# function solve_problem(panels :: Vector{<: Panel2D}, u, num_wake :: Integer, wake_length)
#   wakes     = wake_panels(panels, wake_length, num_wake)
#   φs    = solve_strengths(panels, u, wakes; bound = wake_length)
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