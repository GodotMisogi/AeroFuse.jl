module DoubletSource

export 
solve_strengths, pressure_coefficient, lift_coefficient,
panel_velocities, influence_coefficient

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics

include("../General/math_tools.jl")
# export parts, stencil, midgrad, rotation, inverse_rotation, affine_2D, <<, span, +, *, dot

## Solutions to Laplace's equation
#==========================================================================================#

include("../General/laplace.jl")

export Uniform2D, velocity, source_potential, doublet_potential, grid_data

## Panels
#==========================================================================================#

include("../General/panel.jl")

export Panel2D, split_panels, panel_dist, panel_length, panel_angle, panel_tangent, point1, point2, collocation_point,
make_2Dpanels

#----------------------------Influence coefficient matrices assembly-----------------------------#

doublet_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real) = 
    doublet_potential(strength, affine_2D(x, y, point1(panel)..., panel_angle(panel))..., 0, panel_length(panel))

function source_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real)
    pan = affine_2D(x, y, point1(panel)..., panel_angle(panel))
    source_potential(strength, pan..., 0., panel_length(panel))
end

influence_coefficient(panel_1 :: Panel2D, panel_2 :: Panel2D) = doublet_potential(panel_1, 1., collocation_point(panel_2)...)

doublet_matrix(panels :: Array{Panel2D}) = [ panel_i === panel_j ? 0.5 : influence_coefficient(panel_j, panel_i) for panel_i ∈ panels, panel_j ∈ panels ]

source_matrix(panels :: Array{Panel2D}) = [ source_potential(panel_j, 1., pt...) for pt ∈ collocation_point.(panels), panel_j ∈ panels ]

source_strengths(panels :: Array{Panel2D}, uniform :: Uniform2D) = [ dot(velocity(uniform), normal) for normal ∈ panel_normal.(panels) ]

boundary_condition(panels :: Array{Panel2D}, uniform :: Uniform2D) = - source_matrix(panels) * source_strengths(panels, uniform)

# Morino's velocity Kutta condition
kutta_condition(panels :: Array{<: Panel2D}) = [ 1., -1., zeros(length(panels) - 4)..., 1., -1.]

function wake_vector(panels :: Array{Panel2D}, bound = 1e3)
    lastx, lasty = point2(panels[end])
    woke_panel = Panel2D((lastx, lasty), (bound * lastx, lasty))
    
    [ doublet_potential(woke_panel, 1., pt...) for pt ∈ collocation_point.(panels) ]
end

influence_matrix(panels :: Array{Panel2D}) = [ doublet_matrix(panels)   wake_vector(panels) ;
                                               kutta_condition(panels)'          0.          ]

# vcat(hcat(doublet_matrix(panels), wake_vector(panels)), hcat(kutta_condition(panels)', 0))
                     
boundary_vector(panels :: Array{Panel2D}, uniform :: Uniform2D) = [ boundary_condition(panels, uniform); 0. ]

#---------------------------Solution methods and computations-------------------------------#

panel_velocities(panels :: Array{Panel2D}, uniform :: Uniform2D, doublet_strengths :: Array{<: Any}) = let
    diff_pans = [ panel_dist(panel_1, panel_2) for (panel_1, panel_2) ∈ (collect ∘ eachrow ∘ midgrad)(panels) ]
    diff_strs = [ (str1 - str2) for (str1, str2) ∈ (collect ∘ eachrow ∘ midgrad)(doublet_strengths) ]
    tan_dot_u = [ dot(velocity(uniform), tangent) for tangent ∈ panel_tangent.(panels) ]

    diff_strs ./ diff_pans .+ tan_dot_u 
end

lift_coefficient(cp :: Real, dist_colpoints :: Real, panel_angle :: Real) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength :: Real, speed :: Real) = - 2. * wake_strength / speed

function solve_strengths(panels :: Array{Panel2D}, uniform :: Uniform2D) 
    σs = source_strengths(panels, uniform)
    φs = influence_matrix(panels) \ boundary_vector(panels, uniform)

    φs, σs
end


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