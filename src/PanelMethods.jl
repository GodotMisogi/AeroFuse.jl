module PanelMethods

using LinearAlgebra
using Base.Iterators

⊗(A, B) = kron(A, B)
×(xs, ys) = (collect ∘ zip)(xs' ⊗ (ones ∘ length)(ys), (ones ∘ length)(xs)' ⊗ ys)

# Solutions to Laplace's equation
abstract type Solution end

struct Source2D <: Solution; str :: Float64; x0 :: Float64; y0 :: Float64 end

struct Uniform2D <: Solution; mag :: Float64; ang :: Float64 end

struct Doublet2D <: Solution; str :: Float64; x0 :: Float64; y0 :: Float64 end 

struct Vortex2D <: Solution; str :: Float64; x0 :: Float64; y0 :: Float64 end

# Methods on solutions
velocity(src :: Source2D, x, y) = (src.str / (2π) * (x - src.x0) / ((x - src.x0)^2 + (y - src.y0)^2), str / (2π) * (y - src.y0) / ((x - src.x0)^2 + (y - src.y0)^2))
velocity(uni :: Uniform2D) = let ang = uni.ang * π / 180; (uni.mag * cos(ang), uni.mag * sin(ang)) end
velocity(dub :: Doublet2D, x, y) = (dub.str / (2π) * ((x - dub.x0)^2 - (y - dub.y0)^2) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2, - dub.str / (2π) * 2 * (x - dub.x0) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2)
velocity(vor :: Vortex2D, x, y) = (-vor.str / (2π) * (y - vor.y0) / ((x - vor.x0)^2 + (y - vor.y0)^2), str / (2π) * (x - vor.x0) / ((x - vor.x0)^2 + (y - vor.y0)^2))

potential(src :: Source2D, x, y) = src.str / (4π) * log((x - src.x0)^2 + (y - src.y0)^2)
potential(uni :: Uniform2D, x, y) = let ang = uni.ang * π / 180; uni.mag * (x * cos(ang) + y * sin(ang)) end
potential(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
potential(vor :: Vortex2D, x, y) = vor.str / (2π) * atan(y - vor.y0, x - vor.x0)

stream(src :: Source2D, x, y) = src.str / (2π) * atan(y - src.y0, x - src.x0)
stream(dub :: Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
stream(vor :: Vortex2D, x, y) = -vor.str / (4π) * log((x - vor.x0)^2 + (y - vor.y0)^2)

# Panels for singularity method using Green's third identity
abstract type Panel2D <: Solution end

mutable struct DoubletPanel2D <: Panel2D
    start :: Tuple{Float64, Float64}
    finish :: Tuple{Float64, Float64}
    cp :: Float64
    strength :: Float64
    DoubletPanel2D(start, finish, strength = 0, cp = 0) = new(start, finish, strength, cp)
end

# Methods on panels
collocation_point(panel :: Panel2D) = ((panel.start[1] + panel.finish[1]) / 2, (panel.start[2] + panel.finish[2]) / 2)

function panel_length(panel :: Panel2D) 
    (xs, ys) = panel.start
    (xe, ye) = panel.finish
    return norm([xe - xs, ye - ys])
end

function panel_angle(panel :: Panel2D)
    (xs, ys) = panel.start
    (xe, ye) = panel.finish
    return atan(ye - ys, xe - xs)
end 

panel_tangent(panel :: Panel2D) = rotation(1, 0, -1 * panel_angle(panel))
panel_normal(panel :: Panel2D) = inverse_rotation(0, 1, panel_angle(panel))
panel_location(panel :: Panel2D) = let angle = panel_angle(panel); (π / 2 <= angle <= π) || (-π <= angle <= -π / 2) ? "lower" : "upper" end
split_panels(panels) = [ panel for panel in panels if panel_location(panel) == "upper" ], [ panel for panel in panels if panel_location(panel) == "lower" ]

doublet_potential(xp, yp, len) = -1 / (2π) * (atan(yp, xp - len) - atan(yp, xp - 0))

doublet_velocity(str, x, y, xsl, xel) = str / (2π) .* ( - (y / ((x - xsl)^2 + y^2) - y / ((x - xel)^2 + y^2) ), ( (x - xsl) / ((x - xsl)^2 + y^2) - (x - xel) / ((x - xel)^2 + y^2)) )

function potential(panel :: DoubletPanel2D, x, y)
    len = panel_length(panel)
    angle = panel_angle(panel)
    (xp, yp) = panel_coords(x, y, panel.start..., angle)

    doublet_potential(xp, yp, len)
end

function velocity(panel :: DoubletPanel2D, x0, y0)
    x1, y1 = panel.start
    x2, y2 = panel.finish
    angle = panel_angle(panel)
    x, y = panel_coords(x0, y0, panel.start..., angle)
    u, v = doublet_velocity(panel.strength, x, y, x1, x2)
    
    inverse_rotation(u, v, angle)
end

"""
Computes the distance between two panels.
"""
function dist((p1, p2) :: Tuple{Panel2D,Panel2D}) 
    (xc1, yc1) = collocation_point(p1)
    (xc2, yc2) = collocation_point(p2)
    
    norm([xc2 - xc1, yc2 - yc1])
end

# Solving on list of panels
function solve_strengths(panels :: Array{DoubletPanel2D}, uniform :: Uniform2D) 

    num_panels = length(panels)     # Number of panels
    colpoints = collocation_point.(panels) # Collocation points

    # Doublet submatrix
    doublet_matrix = [ i == j ? 0.5 : potential(panel_j, pt...) for (i, pt) ∈ enumerate(colpoints), (j, panel_j) ∈ enumerate(panels) ]

    # Morino's velocity Kutta condition
    kutta = zeros(num_panels + 1)
    kutta[1] = 1.
    kutta[2] = -1.
    kutta[end-2] = 1.
    kutta[end-1] = -1.
    
    # Wake vector
    (lastx, lasty) = panels[end].finish
    woke_panel = DoubletPanel2D((lastx, lasty), (1e5 * lastx, lasty))
    woke_vector = [ potential(woke_panel, pt...) for pt ∈ colpoints ]

    # Influence matrix with explicit Kutta condition
    influence_matrix = zeros(num_panels + 1, num_panels + 1)
    influence_matrix[1:end-1,1:end-1] = doublet_matrix
    influence_matrix[1:end-1,end] = woke_vector 
    influence_matrix[end,:] = kutta

    # # Boundary condition
    boundary_condition = zeros(num_panels + 1)
    boundary_condition[1:end-1] = [ -potential(uniform, pt...) for pt ∈ colpoints ]

    strengths = influence_matrix \ boundary_condition

    for (panel, strength) in zip(panels, strengths)
        panel.strength = strength
    end

    strengths
end 

function aero_coefficients(panels :: Array{DoubletPanel2D}, uniform :: Uniform2D)
    strengths = solve_strengths(panels, uniform)

    diffpans = dist.(midgrad(panels))
    diffstrs = [ (str2 - str1) for (str2, str1) in (midgrad ∘ take)(strengths, length(panels)) ]

    u = velocity(uniform)
    freestream_speed = norm(u)
    
    tan_vels = diffstrs ./ diffpans

    pressure_coeffs = [ pressure_coefficient((0, vt), freestream_speed) for vt in tan_vels ]

    for (panel, cp) in zip(panels, pressure_coeffs)
        panel.cp = cp
    end
    
    lift_coeff = -sum(pressure_coeffs .* (diffpans ./ 2) .* (cos ∘ panel_angle).(panels))
    kutta_coeff = -2 * last(strengths) / freestream_speed

    lift_coeff, kutta_coeff, pressure_coeffs
end

adjn(xs, n) = (collect ∘ zip)(drop(xs, n), xs)

function parts(xs) 
    adj = adjn(xs, 1)
    
    first(adj), last(adj)
end

function midgrad(xs) 
    (firstpans, lastpans) = parts(xs)
    midpans = adjn(xs, 2)
    
    [firstpans; midpans; lastpans]
end

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at angle_s.
panel_coords(x, y, x_s, y_s, angle_s) = rotation(x - x_s, y - y_s, angle_s)
inverse_rotation(x, y, angle) = (x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle))
rotation(x, y, angle) = (x * cos(angle) + y * sin(angle), -x * sin(angle) + y * cos(angle))

pressure_coefficient(vels, mag) = 1 - (norm(vels))^2 / mag^2

# Performs velocity and potential calculations on a grid
function grid_data(objects :: Array{<:Solution}, xs)
    vels = foldl((v1, v2) -> [ u .+ v for (u, v) in zip(v1, v2) ], [ velocity(object, xs) for object in objects ])
    pots = foldl((v1, v2) -> v1 .+ v2, [ potential(object, xs) for object in objects ])
    
    vels, pots
end

grid_data(object :: Solution, xs) = velocity(object, xs), potential(object, xs)

# Performs velocity and potential computations for an object on a grid
velocity(object :: Solution, xs) = [ velocity(object, x...) for x in xs ] 
potential(object :: Solution, xs) = [ potential(object, x...) for x in xs ]

doublet_panels(coords :: Array{<: Real, 2}) = reverse([ DoubletPanel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ (collect ∘ eachrow)([coords[2:end,:] coords[1:end-1,:]]) ], dims = 1) 

end

# Influence matrix with implicit Kutta condition
# influence_matrix = doublet_matrix
# influence_matrix[:,1] -= woke_vector 
# influence_matrix[:,end] += woke_vector

# Boundary condition
# boundary_condition = [ -potential(uniform, pt...) for pt ∈ colpoints ]

# With sources
# tandotu = [ u.*panel_tangent(panel) for panel in panels ]
# tan_vels = (diffstr ./ diffpans) .+ tandotu

# struct SourcePanel2D <: Panel2D
#     start :: Tuple{Float64, Float64}
#     finish :: Tuple{Float64, Float64} 
# end


# struct VortexPanel2D <: Panel2D
#     start :: Tuple{Float64, Float64}
#     finish :: Tuple{Float64, Float64}
# end

# abstract type Panel3D <: Solution end

# struct SourcePanel3D <: Panel3D
#     start :: Tuple{Float64, Float64, Float64}
#     finish :: Tuple{Float64, Float64 ,Float64}
# end

# struct DoubletPanel3D <: Panel3D
#     start :: Tuple{Float64, Float64, Float64}
#     finish :: Tuple{Float64, Float64 ,Float64}
# end

# struct VortexPanel3D <: Panel3D
#     start :: Tuple{Float64, Float64, Float64}
#     finish :: Tuple{Float64, Float64 ,Float64}
# end