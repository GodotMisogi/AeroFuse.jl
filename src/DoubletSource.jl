module DoubletSource

export 
solve_strengths, pressure_coefficient, lift_coefficient,
panel_velocities, velocity, source_potential, doublet_potential, influence_coefficient,
split_panels, panel_dist, panel_length, panel_angle, panel_tangent, point1, point2, collocation_point,
make_2Dpanels, grid_data, 
Panel2D, Uniform2D, affine_2D

include("MathTools.jl")
include("LaplaceSolutions.jl")

using LinearAlgebra
using Base.Iterators
using StaticArrays
using Statistics
using .MathTools: parts, stencil, midgrad, rotation, inverse_rotation, affine_2D, <<, span, +, *, dot
import Base: +, -, zero

# Crap for automatic differentiation
# zero(:: NTuple{2,<:Real}) = (0., 0.)
# +(::Union{Nothing, Panel2D}, ::Union{Nothing,Panel2D}) = nothing
# zero(:: Nothing) = nothing


#------------------- Solutions to Laplace's equation--------------------#

abstract type Laplace end

struct Uniform2D <: Laplace
    mag :: Real
    ang :: Real 
end

velocity(uni :: Uniform2D) = let ang = deg2rad(uni.ang); (uni.mag * cos(ang), uni.mag * sin(ang)) end

potential(uni :: Uniform2D, x :: Real, y :: Real) = let ang = deg2rad(uni.ang); uni.mag * (x * cos(ang) + y * sin(ang)) end

# Performs velocity and potential calculations on a grid
function grid_data(objects :: Array{<: Laplace}, xs)
    vels = foldl((v1, v2) -> [ u .+ v for (u, v) ∈ zip(v1, v2) ], [ velocity(object, xs) for object ∈ objects ])
    pots = foldl((v1, v2) -> v1 + v2, [ potential(object, xs) for object ∈ objects ])
    
    vels, pots
end

# Performs velocity and potential computations for an object on a grid
grid_data(object :: Laplace, xs) = velocity(object, xs), potential(object, xs)
velocity(object :: Laplace, xs) = map(x -> velocity(object, x...), xs) 
potential(object :: Laplace, xs) = map(x -> potential(object, x...), xs)


abstract type Panel <: Laplace end

struct Panel2D <: Panel
    p1 :: NTuple{2, Real}
    p2 :: NTuple{2, Real}
end

# Methods on panels
point1(p :: Panel) = p.p1
point2(p :: Panel) = p.p2
zero(:: Panel2D) = Panel2D((0., 0.), (0., 0.))


a :: Panel + b :: Panel = Panel2D(point1(a) + point1(b), point2(a) + point2(b))
a :: Panel - b :: Panel = Panel2D(point1(a) - point1(b), point2(a) - point2(b))

collocation_point(panel :: Panel2D) = (point1(panel) .+ point2(panel)) ./ 2
panel_length(panel :: Panel2D) = norm(point2(panel) .- point1(panel))

panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) .- collocation_point(panel_1))
split_panels(panels :: Array{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

make_2Dpanels(coords :: Array{<: Real, 2}) = [ Panel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ eachrow([ coords[2:end,:] coords[1:end-1,:] ]) ][end:-1:1]


function panel_angle(panel :: Panel2D)
    xs, ys = point1(panel) 
    xe, ye = point2(panel)
    
    atan(ye - ys, xe - xs) 
end

panel_tangent(panel :: Panel2D) = rotation(1., 0., -panel_angle(panel))
panel_normal(panel :: Panel2D) = inverse_rotation(0., 1., panel_angle(panel))
panel_location(panel :: Panel2D) = let angle = panel_angle(panel); (π/2 <= angle <= π) || (-π <= angle <= -π/2) ? "lower" : "upper" end

#------Integrated solutions in local panel coordinates for lumped distributions------#

source_potential(str :: Real, x :: Real, z :: Real, x1 :: Real, x2 :: Real) = str / (4π) * ((x - x1) * log((x - x1)^2 + z^2) - (x - x2) * log((x - x2)^2 + z^2) + 2z * (atan(z, x - x2) - atan(z, x - x1)))

source_velocity(str :: Real, x :: Real, z :: Real, x1 :: Real, x2 :: Real) = (str / (4π) * log(((x - x1)^2 + z^2) / ((x - x2)^2 + z^2)), doublet_potential(str, x, z, x1, x2))

doublet_potential(str :: Real, x :: Real, z :: Real, x1 :: Real, x2 :: Real) = str / (2π) * (atan(z, x - x1) - atan(z, x - x2))

doublet_velocity(str :: Real, x :: Real, z :: Real, x1 :: Real, x2 :: Real) = (str / (2π) * - (z / ((x - x1)^2 + z^2) - z / ((x - x2)^2 + z^2) ), str / (2π) * ( (x - x1) / ((x - x1)^2 + z^2) - (x - x2) / ((x - x2)^2 + z^2)))

doublet_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real) = 
    doublet_potential(strength, affine_2D(x, y, point1(panel)..., panel_angle(panel))..., 0, panel_length(panel))

function source_potential(panel :: Panel2D, strength :: Real, x :: Real, y :: Real)
    pan = affine_2D(x, y, point1(panel)..., panel_angle(panel))
    source_potential(strength, pan..., 0., panel_length(panel))
end

#----------------------------Influence coefficient matrices assembly-----------------------------#

influence_coefficient(panel_1 :: Panel2D, panel_2 :: Panel2D) = doublet_potential(panel_1, 1., collocation_point(panel_2)...)

doublet_matrix(panels :: Array{Panel2D}) = [ panel_i === panel_j ? 0.5 : influence_coefficient(panel_j, panel_i) for panel_i ∈ panels, panel_j ∈ panels ]

source_matrix(panels :: Array{Panel2D}) = [ source_potential(panel_j, 1., pt...) for pt ∈ collocation_point.(panels), panel_j ∈ panels ]

source_strengths(panels :: Array{Panel2D}, uniform :: Uniform2D) = [ dot(velocity(uniform), normal) for normal ∈ panel_normal.(panels) ]

boundary_condition(panels :: Array{Panel2D}, uniform :: Uniform2D) = - source_matrix(panels) * source_strengths(panels, uniform)
#  source_strengths(panels, uniform)



# Morino's velocity Kutta condition
kutta_condition(panels :: Array{<: Panel2D}) = vcat([1., -1.], zeros(length(panels) - 4), [1., -1.])

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