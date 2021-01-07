module DoubletSource

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using TimerOutputs

## Non-dimensionalization
#==========================================================================================#

include("../Tools/NonDimensional.jl")
import .NonDimensional: pressure_coefficient

## Math tools
#==========================================================================================#

include("../Tools/MathTools.jl")
using .MathTools: rotation, inverse_rotation, affine_2D, span, midgrad

## Solutions to Laplace's equation
#==========================================================================================#

include("../Tools/Laplace.jl")
import .Laplace: source_potential, doublet_potential, potential

## Panels
#==========================================================================================#

include("../Geometry/PanelGeometry.jl")
using .PanelGeometry: Panel, Panel2D, split_panels, panel_dist, panel_length, panel_normal, panel_angle, panel_normal, panel_tangent, point1, point2, collocation_point, panel_pairs, trans_panel

export Panel, Panel2D, collocation_point

## Matrix assembly
#==========================================================================================#

function doublet_influence(panel_j :: Panel2D, panel_i :: Panel2D)
    xp, yp = trans_panel(panel_j, panel_i)
    doublet_potential(1., xp, yp, 0., panel_length(panel_j))
end

function source_influence(panel_j :: Panel2D, panel_i :: Panel2D)
    xp, yp = trans_panel(panel_j, panel_i)
    source_potential(1., xp, yp, 0., panel_length(panel_j))
end

include("matrix_func.jl")
include("matrix_prealloc.jl")

export solve_problem, solve_strengths, lift_coefficient, solve_case,
panel_velocities, influence_coefficient, doublet_potential, source_potential,
doublet_matrix, wake_vector, source_matrix, boundary_vector, kutta_condition

function solve_strengths(panels :: AbstractVector{<: Panel2D}, u, bound = 1e3) 
    # @timeit "Source Strengths" σs = source_strengths(panels, freestream)
   
    # @timeit "AIC" AIC   = influence_matrix(panels)
    # @timeit "RHS" boco  = boundary_vector(panels, u)


    # Pre-allocated version
    n = length(panels) + 1
    AIC = zeros(n,n)
    boco = zeros(n)

    lastx, lasty = (point2 ∘ last)(panels)
    woke_panel = Panel2D{Float64}((lastx, lasty), (bound * lastx, lasty))

    @timeit "Matrix Assembly (Preallocated)" matrix_assembly!(AIC, boco, panels, woke_panel, u)

    @timeit "Solve AIC" AIC \ boco 
end

function solve_problem(panels :: AbstractVector{<: Panel2D}, u)
    @timeit "Solve System" φs = solve_strengths(panels, u)
        
    @timeit "Lift Coefficient" cl = lift_coefficient(panels, u, φs)

    φs, cl
end

## Force and pressure evaluations
#==========================================================================================#

function panel_velocities(panels :: AbstractVector{<: Panel2D}, u, doublet_strengths :: AbstractVector{<: Real})
    @timeit "Panel Pairs" diff_pans = panel_pairs(panels)
    @timeit "Strength Diffs" diff_strs = -diff(midgrad(doublet_strengths), dims = 2)
    @timeit "Tangential Velocities" tan_dot_u = dot.(Ref(u), panel_tangent.(panels))

    @timeit "Sum Velocities" @. diff_strs / diff_pans + tan_dot_u 
end

lift_coefficient(cp :: Real, dist_colpoints :: Real, panel_angle :: Real) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength :: Real, speed :: Real) = - 2. * wake_strength / speed

function pressure_coefficient(panels :: AbstractVector{<: Panel2D}, φs :: AbstractVector{<: Real}, u) 
    @timeit "Panel Velocities" panel_vels = panel_velocities(panels, u, φs[1:end-1])
    pressure_coefficient.(norm(u), panel_vels)
end

function lift_coefficient(panels :: AbstractVector{<: Panel2D}, u, φs :: AbstractVector{<: Real}) 
    @timeit "Pressure Coefficients" cps = pressure_coefficient(panels, φs, u)
    @timeit "Lift Coefficients" cls = lift_coefficient.(cps, panel_pairs(panels) / 2, panel_angle.(panels))
    cl = sum(cls)
end

# Make panels for plotting
# dub_src_panels = DoubletSourcePanel2D.(panels, φs[1:end-1], σs, pressure_coeffs)


#-----------------------Doublet-source panel setup for plotting---------------------------#

# struct DoubletSourcePanel2D <: Panel
#     coords :: Panel2D
#     doublet_strength :: Float64
#     source_strength :: Float64
#     cp :: Float64
# end

# split_panels(panels :: Array{DoubletSourcePanel2D}) = collect.(span(panel -> panel_location(panel.coords) == "upper", panels))

    
# potential(panel :: DoubletSourcePanel2D, x :: Real, y :: Real) = let len = panel_length(panel.coords);
#     doublet_potential(panel.doublet_strength, x, y, 0, len) + source_potential(panel.source_strength, x, y, 0, len)
# end

# function velocity(panel :: DoubletSourcePanel2D, x :: Real, y :: Real)
#     x1, y1 = point1(panel.coords)
#     x2, y2 = point2(panel.coords)
#     α = panel_angle(panel.coords)
#     xp, yp = affine_2D(x, y, point1(panel.coords)..., α)
#     u, v = doublet_velocity(panel.doublet_strength, xp, yp, x1, x2) .+ source_velocity(panel.source_strength, xp, yp, x1, x2)
    
#     inverse_rotation(u, v, α)
# end

end