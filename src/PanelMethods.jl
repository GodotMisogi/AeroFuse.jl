module PanelMethods

using LinearAlgebra
using Base.Iterators

⊗(A, B) = kron(A, B)
×(xs, ys) = (collect ∘ zip)(xs' ⊗ (ones ∘ length)(ys), (ones ∘ length)(xs)' ⊗ ys)

span(pred, iter) = (takewhile(pred, iter), dropwhile(pred, iter))
lisa(pred, iter) = span(!pred, iter)

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
    strength :: Float64
    cp :: Float64
    DoubletPanel2D(start, finish, strength = 1, cp = 0) = new(start, finish, strength, cp)
end

doublet_panels(coords :: Array{<: Real, 2}) = reverse([ DoubletPanel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ (collect ∘ eachrow)([coords[2:end,:] coords[1:end-1,:]]) ], dims = 1) 

mutable struct DoubletSourcePanel2D <: Panel2D
    start :: Tuple{Float64, Float64}
    finish :: Tuple{Float64, Float64}
    doublet_strength :: Float64
    source_strength :: Float64
    cp :: Float64
    DoubletSourcePanel2D(start, finish, doublet_strength = 1, source_strength = 1, cp = 0) = new(start, finish, doublet_strength, source_strength, cp)
end

doublet_source_panels(coords :: Array{<: Real, 2}) = reverse([ DoubletSourcePanel2D((xs, ys), (xe, ye)) for (xs, ys, xe, ye) ∈ (collect ∘ eachrow)([coords[2:end,:] coords[1:end-1,:]]) ], dims = 1) 

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
split_panels(panels) = span(panel -> panel_location(panel) == "lower", panels)

"""
Computes the distance between two panels.
"""
function dist((p1, p2) :: Tuple{Panel2D,Panel2D}) 
    (xc1, yc1) = collocation_point(p1)
    (xc2, yc2) = collocation_point(p2)
    
    norm([xc2 - xc1, yc2 - yc1])
end

source_potential(str, x, z, x1, x2) = str / (4π) * ((x - x1) * log((x - x1)^2 + z^2) - (x - x2) * log((x - x2)^2 + z^2) + 2z * (atan(z, x - x2) - atan(z, x - x1)))

source_velocity(str, x, z, x1, x2) = ( str / (4π) * log(((x - x1)^2 + z^2) / ((x - x2)^2 + z^2)), doublet_potential(str, x, z, x1, x2) )

doublet_potential(str, x, z, x1, x2) = str / (2π) * (atan(z, x - x1) - atan(z, x - x2))

doublet_velocity(str, x, z, x1, x2) = str / (2π) .* ( - (z / ((x - x1)^2 + z^2) - z / ((x - x2)^2 + z^2) ), ( (x - x1) / ((x - x1)^2 + z^2) - (x - x2) / ((x - x2)^2 + z^2)) )

function influence_matrix(panels)
    
    num_panels = length(panels)

    # Influence matrix with explicit Kutta condition
    influence_matrix = zeros(num_panels + 1, num_panels + 1)
    influence_matrix[1:end-1,1:end-1] = doublet_matrix(panels)
    influence_matrix[1:end-1,end] = wake_vector(panels)
    influence_matrix[end,:] = morino_condition(panels)

    influence_matrix
end

function boundary_vector(panels, uniform)
    # Boundary condition
    boundary_vector = zeros(length(panels) + 1)
    boundary_vector[1:end-1] = boundary_condition(panels, uniform)

    boundary_vector
end


#---------------DOUBLET PANEL METHOD-----------------#


function potential(panel :: DoubletPanel2D, x, y)
    len = panel_length(panel)
    angle = panel_angle(panel)
    xp, yp = panel_coords(x, y, panel.start..., angle)

    doublet_potential(panel.strength, xp, yp, 0, len)
end

function velocity(panel :: DoubletPanel2D, x, y)
    x1, y1 = panel.start
    x2, y2 = panel.finish
    angle = panel_angle(panel)
    xp, yp = panel_coords(x, y, panel.start..., angle)
    u, v = doublet_velocity(panel.strength, xp, yp, x1, x2)
    
    inverse_rotation(u, v, angle)
end

doublet_matrix(panels :: Array{DoubletPanel2D}) = [ i == j ? 0.5 : potential(panel_j, pt...) for (i, pt) ∈ enumerate(collocation_point.(panels)), (j, panel_j) ∈ enumerate(panels) ]

boundary_condition(panels :: Array{DoubletPanel2D}, uniform) = [ -potential(uniform, pt...) for pt ∈ collocation_point.(panels) ]

function wake_vector(panels :: Array{DoubletPanel2D})
    # Wake vector
    lastx, lasty = panels[end].finish
    woke_panel = DoubletPanel2D(panels[end].finish, (1e5 * lastx, lasty))
    
    [ potential(woke_panel, pt...) for pt ∈ collocation_point.(panels) ]
end 

function solve_strengths(panels :: Array{DoubletPanel2D}, uniform :: Uniform2D) 

    strengths = influence_matrix(panels) \ boundary_vector(panels, uniform)

    for (panel, strength) in zip(panels, strengths)
        panel.strength = strength
    end

    strengths
end

function panel_velocities(panels :: Array{DoubletPanel2D})
    strengths = [ panel.strength for panel in panels ]
    diff_pans = dist.(midgrad(panels))
    diff_strs = [ (str2 - str1) for (str2, str1) in (midgrad ∘ take)(strengths, length(panels)) ]

    panel_vels = diff_strs ./ diff_pans
end

lift_coefficient(cp, dist_colpoints, panel_angle) = - cp * dist_colpoints * cos(panel_angle)

lift_coefficient(wake_strength, speed) = - 2 * wake_strength / speed

function aero_coefficients(panels :: Array{DoubletPanel2D}, uniform :: Uniform2D)
    strengths = solve_strengths(panels, uniform)

    freestream_speed = (norm ∘ velocity)(uniform)

    pressure_coeffs = [ pressure_coefficient((0, vt), freestream_speed) for vt in panel_velocities(panels) ]

    lift_coeff = sum(lift_coefficient.(pressure_coeffs, dist.(midgrad(panels)) ./ 2, panel_angle.(panels)))

    kutta_coeff = lift_coefficient(last(strengths), freestream_speed)

    for (panel, cp) in zip(panels, pressure_coeffs)
        panel.cp = cp
    end

    lift_coeff, kutta_coeff, pressure_coeffs
end

#---------------DOUBLET-SOURCE PANEL METHOD-----------------#

function doublet_potential(panel :: DoubletSourcePanel2D, x, y)
    len = panel_length(panel)
    angle = panel_angle(panel)
    (xp, yp) = panel_coords(x, y, panel.start..., angle)

    doublet_potential(panel.doublet_strength, xp, yp, 0, len)
end

function source_potential(panel :: DoubletSourcePanel2D, x, y)
    len = panel_length(panel)
    angle = panel_angle(panel)
    (xp, yp) = panel_coords(x, y, panel.start..., angle)

    source_potential(panel.source_strength, xp, yp, 0, len)
end

function potential(panel :: DoubletSourcePanel2D, x, y) 
    len = panel_length(panel)
    angle = panel_angle(panel)
    xp, yp = panel_coords(x, y, panel.start..., angle)

    doublet_potential(panel, x, y, 0, len) + souce_potential(panel, x, y, 0, len)
end

function velocity(panel :: DoubletSourcePanel2D, x, y)
    x1, y1 = panel.start
    x2, y2 = panel.finish
    angle = panel_angle(panel)
    xp, yp = panel_coords(x, y, panel.start..., angle)
    u, v = doublet_velocity(panel.doublet_strength, xp, yp, x1, x2) .+ source_velocity(panel.source_strength, xp, yp, x1, x2)
    
    inverse_rotation(u, v, angle)
end

doublet_matrix(panels :: Array{DoubletSourcePanel2D}) = [ i == j ? 0.5 : doublet_potential(panel_j, pt...) for (i, pt) ∈ enumerate(collocation_point.(panels)), (j, panel_j) ∈ enumerate(panels) ]

source_matrix(panels :: Array{DoubletSourcePanel2D}) = [ source_potential(panel_j, pt...) for pt ∈ collocation_point.(panels), panel_j ∈ panels ]
source_strengths(panels :: Array{DoubletSourcePanel2D}, uniform :: Uniform2D) = [ sum(normal .* velocity(uniform)) for normal in panel_normal.(panels) ]

boundary_condition(panels :: Array{DoubletSourcePanel2D}, uniform :: Uniform2D) = - source_matrix(panels) * source_strengths(panels, uniform)

function wake_vector(panels :: Array{DoubletSourcePanel2D})
    # Wake vector
    lastx, lasty = panels[end].finish
    woke_panel = DoubletSourcePanel2D(panels[end].finish, (1e5 * lastx, lasty))
    
    [ doublet_potential(woke_panel, pt...) for pt ∈ collocation_point.(panels) ]
end

# Morino's velocity Kutta condition
function morino_condition(panels :: Array{<: Panel2D})
    num_panels = length(panels)
    kutta = zeros(num_panels + 1)
    kutta[1] = 1.
    kutta[2] = -1.
    kutta[end-2] = 1.
    kutta[end-1] = -1.

    kutta
end

# Solving on list of panels
function panel_velocities(panels :: Array{DoubletSourcePanel2D}, uniform)
    strengths = [ panel.doublet_strength for panel in panels ]
    diff_pans = dist.(midgrad(panels))
    diff_strs = [ (str2 - str1) for (str2, str1) in (midgrad ∘ take)(strengths, length(panels)) ]
    tan_dot_u = [ sum(velocity(uniform) .* tangent) for tangent in panel_tangent.(panels) ]

    panel_vels = diff_strs ./ diff_pans .+ tan_dot_u 
end

function solve_strengths(panels :: Array{DoubletSourcePanel2D}, uniform :: Uniform2D) 

    dub_strengths = influence_matrix(panels) \ boundary_vector(panels, uniform)
    src_strengths = source_strengths(panels, uniform)

    for (panel, doublet_strength, source_strength) in zip(panels, dub_strengths, src_strengths)
        panel.doublet_strength = doublet_strength
        panel.source_strength = source_strength
    end

    dub_strengths, src_strengths
end 

function aero_coefficients(panels :: Array{DoubletSourcePanel2D}, uniform :: Uniform2D)
    dub_strengths, src_strengths = solve_strengths(panels, uniform)

    freestream_speed = (norm ∘ velocity)(uniform)
    
    pressure_coeffs = [ pressure_coefficient((0, vt), freestream_speed) for vt in panel_velocities(panels, uniform) ]

    lift_coeff = sum(lift_coefficient.(pressure_coeffs, dist.(midgrad(panels)) ./ 2, panel_angle.(panels)))

    for (panel, cp) in zip(panels, pressure_coeffs)
        panel.cp = cp
    end

    lift_coeff, pressure_coeffs
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