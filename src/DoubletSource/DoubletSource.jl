module DoubletSource


using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using TimerOutputs

include("../Tools/MathTools.jl")
using .MathTools: rotation, inverse_rotation, affine_2D, span, midgrad

## Solutions to Laplace's equation
#==========================================================================================#

include("../Tools/laplace.jl")

export Uniform2D, velocity, source_potential, doublet_potential, grid_data

## Panels
#==========================================================================================#

include("../Geometry/PanelGeometry.jl")
using .PanelGeometry: Panel, Panel2D, split_panels, panel_dist, panel_length, panel_normal, panel_angle, panel_normal, panel_tangent, point1, point2, collocation_point, panel_pairs

export Panel, Panel2D, collocation_point

## Matrix assembly
#==========================================================================================#

export solve_strengths, pressure_coefficient, lift_coefficient, solve_case,
panel_velocities, influence_coefficient, doublet_potential, source_potential,
doublet_matrix, wake_vector, source_matrix, boundary_vector, kutta_condition

doublet_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real) = 
    @timeit "Doublet Potential" doublet_potential(strength, affine_2D(x, y, point1(panel)..., panel_angle(panel))..., 0, panel_length(panel))

function source_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real)
    pan = affine_2D(x, y, point1(panel)..., panel_angle(panel))
    @timeit "Source Potential" source_potential(strength, pan..., 0., panel_length(panel))
end

doublet_influence(panel_1 :: Panel2D, panel_2 :: Panel2D) = doublet_potential(panel_1, 1., collocation_point(panel_2)...)

source_influence(panel_1 :: Panel2D, panel_2 :: Panel2D) = source_potential(panel_1, 1., collocation_point(panel_2)...)

"""
    doublet_matrix(panels_1, panels_2)

Computes the matrix of doublet potential influence coefficients between pairs of panels_1 and panels_2.
"""
doublet_matrix(panels_1 :: AbstractVector{Panel2D}, panels_2 :: AbstractVector{Panel2D}) = 
    map(pans -> pans[1] === pans[2] ? 0.5 : doublet_influence(pans...), product(panels_1, panels_2))'

"""
    source_matrix(panels_1, panels_2)

Computes the matrix of source potential influence coefficients between pairs of panels_1 and panels_2.
"""
source_matrix(panels_1 :: AbstractVector{Panel2D}, panels_2 :: AbstractVector{Panel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

source_strengths(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D) = [ dot(velocity(freestream), normal) for normal ∈ panel_normal.(panels) ]

boundary_condition(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D) = - source_matrix(panels, panels) * source_strengths(panels, freestream)

# Morino's velocity Kutta condition
kutta_condition(panels :: AbstractVector{Panel2D}) = [ 1., -1., zeros(length(panels) - 4)..., 1., -1.]

function wake_vector(panels :: AbstractVector{Panel2D}, bound = 1e3)
    lastx, lasty = point2(panels[end])
    woke_panel = Panel2D((lastx, lasty), (bound * lastx, lasty))
    
    [ doublet_potential(woke_panel, 1., pt...) for pt ∈ collocation_point.(panels) ]
end

function influence_matrix(panels :: AbstractVector{Panel2D}) 
    @timeit "Doublet Matrix" dub_mat = doublet_matrix(panels, panels)
    @timeit "Wake Vector" wake_vec = wake_vector(panels)
    @timeit "Kutta Condition" kutta = kutta_condition(panels)'
    @timeit "Matrix Assembly" mat = [ dub_mat   wake_vec;
                                      kutta        0.   ]
end
              
boundary_vector(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D) = [ boundary_condition(panels, freestream); 0. ]

function solve_strengths(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D) 
    # @timeit "Source Strengths" σs = source_strengths(panels, freestream)
    @timeit "AIC" AIC = influence_matrix(panels)
    @timeit "RHS" boco = boundary_vector(panels, freestream)
    @timeit "Solve AIC" φs = AIC \ boco 

    φs
end

## Force and pressure evaluations
#==========================================================================================#

function panel_velocities(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D, doublet_strengths :: AbstractVector{<: Real})
    @timeit "Panel Pairs" diff_pans = panel_pairs(panels)
    @timeit "Strength Diffs" diff_strs = -diff(midgrad(doublet_strengths), dims = 2)
    @timeit "Tangential Velocities" tan_dot_u = [ dot(velocity(freestream), tangent) for tangent ∈ panel_tangent.(panels) ]

    @timeit "Sum Velocities" diff_strs ./ diff_pans .+ tan_dot_u 
end

lift_coefficient(cp :: Real, dist_colpoints :: Real, panel_angle :: Real) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength :: Real, speed :: Real) = - 2. * wake_strength / speed

"""
Computes the incompressible pressure coefficient given a magnitude and a velocity vector.
"""
pressure_coefficient(mag, vels) = 1 - norm(vels)^2 / mag^2

function pressure_coefficient(panels :: AbstractVector{Panel2D}, φs :: AbstractVector{<: Real}, freestream :: Uniform2D) 
    speed = (norm ∘ velocity)(freestream)
    @timeit "Panel Velocities" panel_vels = panel_velocities(panels, freestream, φs[1:end-1])
    pressure_coefficient.(speed, panel_vels)
end

function lift_coefficient(panels :: AbstractVector{Panel2D}, freestream :: Uniform2D, φs :: AbstractVector{<: Real}) 
    @timeit "Pressure Coefficients" cps = pressure_coefficient(panels, φs, freestream)
    @timeit "Lift Coefficients" cls = lift_coefficient.(cps, panel_pairs(panels) ./ 2, panel_angle.(panels))
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