module DoubletSource

export solve_case, pressure_coefficient, make_panels, Uniform2D, grid_data, split_panels, collocation_point

include("MathTools.jl")
include("LaplaceSolutions.jl")
using LinearAlgebra
using Base.Iterators
using StaticArrays
using .MathTools: parts, stencil, midgrad, rotation, inverse_rotation, affine_2D, <<, span, dot, structtolist

pressure_coefficient(mag, vels) = 1 - norm(vels)^2 / mag^2


#------------------- Solutions to Laplace's equation--------------------#

abstract type Laplace end

struct Source2D <: Laplace
    str :: Float64
    x0 :: Float64
    y0 :: Float64 
end

struct Uniform2D <: Laplace
    mag :: Float64
    ang :: Float64 
end

struct Doublet2D <: Laplace
    str :: Float64
    x0 :: Float64
    y0 :: Float64 
end 

struct Vortex2D <: Laplace
    str :: Float64
    x0 :: Float64
    y0 :: Float64 
end

velocity(src :: Source2D, x, y) = (src.str / (2π) * (x - src.x0) / ((x - src.x0)^2 + (y - src.y0)^2), str / (2π) * (y - src.y0) / ((x - src.x0)^2 + (y - src.y0)^2))
velocity(uni :: Uniform2D) = let ang = deg2rad(uni.ang); (uni.mag * cos(ang), uni.mag * sin(ang)) end
velocity(dub :: Doublet2D, x, y) = (dub.str / (2π) * ((x - dub.x0)^2 - (y - dub.y0)^2) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2, - dub.str / (2π) * 2 * (x - dub.x0) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2)
velocity(vor :: Vortex2D, x, y) = (-vor.str / (2π) * (y - vor.y0) / ((x - vor.x0)^2 + (y - vor.y0)^2), str / (2π) * (x - vor.x0) / ((x - vor.x0)^2 + (y - vor.y0)^2))

potential(src :: Source2D, x, y) = src.str / (4π) * log((x - src.x0)^2 + (y - src.y0)^2)
potential(uni :: Uniform2D, x, y) = let ang = deg2rad(uni.ang); uni.mag * (x * cos(ang) + y * sin(ang)) end
potential(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
potential(vor :: Vortex2D, x, y) = vor.str / (2π) * atan(y - vor.y0, x - vor.x0)

stream(src :: Source2D, x, y) = src.str / (2π) * atan(y - src.y0, x - src.x0)
stream(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
stream(vor :: Vortex2D, x, y) = -vor.str / (4π) * log((x - vor.x0)^2 + (y - vor.y0)^2)

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

#--------------------------Vector spaces------------------------------#

abstract type Point end

struct Point2D <: Point
    x :: Real; y :: Real;
end

struct Point3D <: Point
    x :: Real; y :: Real; z :: Real;
end

abstract type Panel <: Laplace end


struct Panel2D <: Panel
    p1 :: Tuple{Float64, Float64}
    p2 :: Tuple{Float64, Float64}
end

# Methods on panels in N dimensions
collocation_point(panel :: Panel2D) = (panel.p1 .+ panel.p2) ./ 2
panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) .- collocation_point(panel_1))
split_panels(panels :: Array{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

make_panels(coords :: Array{<: Real, 2}) = [ Panel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ eachrow([ coords[2:end,:] coords[1:end-1,:] ]) ][end:-1:1]

panel_length(panel :: Panel2D) = norm(panel.p2 .- panel.p1)
panel_angle(panel :: Panel2D) = let (xs, ys) = panel.p1, (xe, ye) = panel.p2; atan(ye - ys, xe - xs) end
panel_tangent(panel :: Panel2D) = rotation(1, 0, -1 * panel_angle(panel))
panel_normal(panel :: Panel2D) = inverse_rotation(0, 1, panel_angle(panel))
panel_location(panel :: Panel2D) = let angle = panel_angle(panel); (π/2 <= angle <= π) || (-π <= angle <= -π/2) ? "lower" : "upper" end

struct Panel3D <: Panel
    p1 :: Tuple{Float64, Float64, Float64}
    p2 :: Tuple{Float64, Float64, Float64}
    p3 :: Tuple{Float64, Float64, Float64}
    p4 :: Tuple{Float64, Float64, Float64}
end

collocation_point(panel :: Panel3D) = (panel.p1 .+ panel.p2 .+ panel.p3 .+ panel.p4) ./ 4
panel_normal(panel :: Panel3D) = cross(panel.p2 .- panel.p1, panel.p3 .- panel.p2)

struct Line
    r1 :: Point
    r2 :: Point
end

#-----------------------Doublet-source panel method---------------------------#

struct DoubletSourcePanel2D <: Panel
    coords :: Panel2D
    doublet_strength :: Float64
    source_strength :: Float64
    cp :: Float64
end

split_panels(panels :: Array{DoubletSourcePanel2D}) = collect.(span(panel -> panel_location(panel.coords) == "upper", panels))

#------Integrated solutions in local panel coordinates for lumped distributions------#

source_potential(str, x, z, x1, x2) = str / (4π) * ((x - x1) * log((x - x1)^2 + z^2) - (x - x2) * log((x - x2)^2 + z^2) + 2z * (atan(z, x - x2) - atan(z, x - x1)))

source_velocity(str, x, z, x1, x2) = ( str / (4π) * log(((x - x1)^2 + z^2) / ((x - x2)^2 + z^2)), doublet_potential(str, x, z, x1, x2) )

doublet_potential(str, x, z, x1, x2) = str / (2π) * (atan(z, x - x1) - atan(z, x - x2))

doublet_velocity(str, x, z, x1, x2) = str / (2π) .* ( - (z / ((x - x1)^2 + z^2) - z / ((x - x2)^2 + z^2) ), ( (x - x1) / ((x - x1)^2 + z^2) - (x - x2) / ((x - x2)^2 + z^2)) )

doublet_potential(panel :: Panel2D, strength, x, y) = 
    doublet_potential(strength, affine_2D(x, y, panel.p1..., panel_angle(panel))..., 0, panel_length(panel))

source_potential(panel :: Panel2D, strength, x, y) =
    source_potential(strength, affine_2D(x, y, panel.p1..., panel_angle(panel))..., 0, panel_length(panel))

    
potential(panel :: DoubletSourcePanel2D, x, y) = let len = panel_length(panel.coords);
    doublet_potential(panel.doublet_strength, x, y, 0, len) + source_potential(panel.source_strength, x, y, 0, len)
end

function velocity(panel :: DoubletSourcePanel2D, x, y)
    x1, y1 = panel.coords.p1
    x2, y2 = panel.coords.p2
    α = panel_angle(panel.coords)
    xp, yp = affine_2D(x, y, panel.coords.p1..., α)
    u, v = doublet_velocity(panel.doublet_strength, xp, yp, x1, x2) .+ source_velocity(panel.source_strength, xp, yp, x1, x2)
    
    inverse_rotation(u, v, α)
end

#----------------------------Influence coefficient matrices assembly-----------------------------#

influence_coefficient(panel_1 :: Panel2D, panel_2 :: Panel2D) = doublet_potential(panel_1, 1, collocation_point(panel_2)...)

doublet_matrix(panels :: Array{Panel2D}) = [ panel_i === panel_j ? 0.5 : influence_coefficient(panel_j, panel_i) for panel_i ∈ panels, panel_j ∈ panels ]

source_matrix(panels :: Array{Panel2D}) = [ source_potential(panel_j, 1, pt...) for pt ∈ collocation_point.(panels), panel_j ∈ panels ]

source_strengths(panels :: Array{Panel2D}, uniform :: Uniform2D) = [ dot(velocity(uniform), normal) for normal ∈ panel_normal.(panels) ]

boundary_condition(panels :: Array{Panel2D}, uniform :: Uniform2D) = - source_matrix(panels) * source_strengths(panels, uniform)

# Morino's velocity Kutta condition
kutta_condition(panels :: Array{<: Panel2D}) = [ 1, -1, zeros(length(panels) - 4)..., 1, -1 ]

function wake_vector(panels :: Array{Panel2D}, bound = 1e5)
    lastx, lasty = panels[end].p2
    woke_panel = Panel2D(panels[end].p2, (bound * lastx, lasty))
    
    [ doublet_potential(woke_panel, 1, pt...) for pt ∈ collocation_point.(panels) ]
end

influence_matrix(panels :: Array{Panel2D}) = 
    [ doublet_matrix(panels)   wake_vector(panels) ;
      kutta_condition(panels)'          0          ]

boundary_vector(panels :: Array{Panel2D}, uniform :: Uniform2D) = [ boundary_condition(panels, uniform); 0 ]

#---------------------------Solution methods and computations-------------------------------#

function panel_velocities(panels :: Array{Panel2D}, uniform :: Uniform2D, doublet_strengths)
    diff_pans = [ panel_dist(panel_1, panel_2) for (panel_1, panel_2) ∈ (eachrow ∘ midgrad)(panels) ]
    diff_strs = [ (str1 - str2) for (str1, str2) ∈ (eachrow ∘ midgrad)(doublet_strengths) ]
    tan_dot_u = [ dot(velocity(uniform), tangent) for tangent ∈ panel_tangent.(panels) ]

    diff_strs ./ diff_pans .+ tan_dot_u 
end

lift_coefficient(cp, dist_colpoints, panel_angle) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength, speed) = - 2 * wake_strength / speed

function solve_strengths(panels :: Array{Panel2D}, uniform :: Uniform2D) 
    src_strengths = source_strengths(panels, uniform)
    dub_strengths = influence_matrix(panels) \ boundary_vector(panels, uniform)

    dub_strengths, src_strengths
end

function solve_case(panels :: Array{Panel2D}, uniform :: Uniform2D)

    # Compute doublet and source strengths
    dub_strengths, src_strengths = solve_strengths(panels, uniform)
    freestream_speed = (norm ∘ velocity)(uniform)
    pressure_coeffs = pressure_coefficient.(freestream_speed, panel_velocities(panels, uniform, dub_strengths[1:end-1]))

    # Make panels for plotting
    dub_src_panels = DoubletSourcePanel2D.(panels, dub_strengths[1:end-1], src_strengths, pressure_coeffs)

    # Compute lift coefficient
    diff_pans = [ panel_dist(panel_1, panel_2) for (panel_1, panel_2) ∈ (eachrow ∘ midgrad)(panels) ]
    cl = sum(lift_coefficient.(pressure_coeffs, diff_pans ./ 2, panel_angle.(panels)))

    dub_src_panels, cl
end

#------------------------------------------------------------------#

end