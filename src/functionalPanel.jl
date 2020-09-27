module AeroModules

using LinearAlgebra
using Base.Iterators
using Interpolations

# Solutions to Laplace's equation
abstract type Solution end

struct Source2D <: Solution; str::Float64; x0::Float64; y0::Float64 end

struct Uniform2D <: Solution; mag::Float64; ang::Float64 end

struct Doublet2D <: Solution; str::Float64; x0::Float64; y0::Float64 end 

struct Vortex2D <: Solution; str::Float64; x0::Float64; y0::Float64 end

# Methods on solutions
velocity(src::Source2D, x, y) = (src.str / (2π) * (x - src.x0) / ((x - src.x0)^2 + (y - src.y0)^2), str / (2π) * (y - src.y0) / ((x - src.x0)^2 + (y - src.y0)^2))
velocity(uni::Uniform2D) = let ang = uni.ang * π / 180; (uni.mag * cos(ang), uni.mag * sin(ang)) end
velocity(dub::Doublet2D, x, y) = (-dub.str / (2π) * ((x - dub.x0)^2 - (y - dub.y0)^2) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2, -str / (2π) * 2 * (x - dub.x0) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)^2)
velocity(vor::Vortex2D, x, y) = (-vor.str / (2π) * (y - vor.y0) / ((x - vor.x0)^2 + (y - vor.y0)^2), str / (2π) * (x - vor.x0) / ((x - vor.x0)^2 + (y - vor.y0)^2))

potential(src::Source2D, x, y) = src.str / (4π) * log((x - src.x0)^2 + (y - src.y0)^2)
potential(uni::Uniform2D, x, y) = let ang = uni.ang * π / 180; uni.mag * (x * cos(ang) + y * sin(ang)) end
potential(dub::Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
potential(vor::Vortex2D, x, y) = vor.str / (2π) * atan(y - vor.y0, x - vor.x0)

stream(src::Source2D, x, y) = src.str / (2π) * atan(y - src.y0, x - src.x0)
stream(dub::Doublet2D, x, y) = -dub.str / (2π) * (y - dub.y0) / ((x - dub.x0)^2 + (y - dub.y0)^2)
stream(vor::Vortex2D, x, y) = -vor.str / (4π) * log((x - vor.x0)^2 + (y - vor.y0)^2)

# Panels for singularity method using Green's third identity
abstract type Panel2D <: Solution end

struct SourcePanel2D <: Panel2D; start::Tuple{Float64,Float64}; finish::Tuple{Float64,Float64} end

struct DoubletPanel2D <: Panel2D; start::Tuple{Float64,Float64}; finish::Tuple{Float64,Float64} end

struct VortexPanel2D <: Panel2D; start::Tuple{Float64,Float64}; finish::Tuple{Float64,Float64} end

abstract type Panel3D <: Solution end

struct SourcePanel3D <: Panel3D; start::Tuple{Float64,Float64,Float64}; finish::Tuple{Float64,Float64,Float64} end

struct DoubletPanel3D <: Panel3D; start::Tuple{Float64,Float64,Float64}; finish::Tuple{Float64,Float64,Float64} end

struct VortexPanel3D <: Panel3D; start::Tuple{Float64,Float64,Float64}; finish::Tuple{Float64,Float64,Float64} end

# Methods on panels
collocation_point(panel::Panel2D) = ((panel.start[1] + panel.finish[1]) / 2, (panel.start[2] + panel.finish[2]) / 2)

function panel_length(panel::Panel2D) 
    (xs, ys) = panel.start
    (xe, ye) = panel.finish
    return norm([xe - xs, ye - ys])
end

function panel_angle(panel::Panel2D)
    (xs, ys) = panel.start
    (xe, ye) = panel.finish
    return atan(ye - ys, xe - xs)
end 

panel_tangent(panel::Panel2D) = rotation(1, 0, -1 * panel_angle(panel))
panel_normal(panel::Panel2D) = inverse_rotation(0, 1, panel_angle(panel))
panel_location(panel::Panel2D) = let angle = panel_angle(panel); (π / 2 <= angle <= π) || (-π <= angle <= -π / 2) ? "lower" : "upper" end

function influence_potential(panel::DoubletPanel2D, x, y)
    # Analytical solution
    (x0, y0) = panel.start
    len = panel_length(panel)
    angle = panel_angle(panel)
    (xp, yp) = panelCoords(x, y, x0, y0, angle)

    return -1 / (2π) * (atan(yp, xp - len) - atan(yp, xp - 0))
end

function dist2((p1, p2)::Tuple{Panel2D,Panel2D}) 
    (xc1, yc1) = collocation_point(p1)
    (xc2, yc2) = collocation_point(p2)
    return norm([ xc2 - xc1, yc2 - yc1 ])
end

# Solving on list of panels
function solve_strengths(panels::Array{DoubletPanel2D,1}, uniform::Uniform2D, sources = false) 

    num_panels = length(panels)     # Number of panels

    # Doublet matrix
    doubletMatrix = [ i == j ? 0.5 : influence_potential(panel_j, xc, yc) for (i, (xc, yc)) in enumerate(collocation_point.(panels)), (j, panel_j) in enumerate(panels) ]

    # Morino's velocity Kutta condition
    kutta = zeros(num_panels + 1)
    kutta[1] = 1.
    kutta[2] = -1.
    kutta[end - 2] = 1.
    kutta[end - 1] = -1.
    
    # Wake vector
    (lastx, lasty) = panels[end].finish
    wokePanel = DoubletPanel2D((lastx, lasty), (100000 * lastx, lasty))
    wokeVector = [ influence_potential(wokePanel, xc, yc) for (xc, yc) in collocation_point.(panels) ]

    # Influence matrix with explicit Kutta condition
    influenceMatrix = zeros(num_panels + 1, num_panels + 1)
    influenceMatrix[1:end - 1, 1:end - 1] = doubletMatrix
    influenceMatrix[1:end - 1, end] = wokeVector 
    influenceMatrix[end, :] = kutta

    # Boundary condition
    u = velocity(uniform)
    boundaryCondition = zeros(num_panels + 1)
    boundaryCondition[1:end - 1] = [ -potential(uniform, xc, yc) for (xc, yc) in collocation_point.(panels) ]
    
    return influenceMatrix \ boundaryCondition
end 

function aerocoefficients(panels::Array{DoubletPanel2D,1}, uniform::Uniform2D)
    strengths = solve_strengths(panels, uniform)

    diffpans = dist2.(midgrad(panels))
    diffstrs = [ (str2 - str1) for (str2, str1) in midgrad(take(strengths, length(panels))) ]

    u = velocity(uniform)
    margaret = norm(u)
    
    # With sources
    # tandotu = [ u.*panelTangent(panel) for panel in panels ]
    # tanVels = map((ds, dp, ut) -> ds/dp + ut, zip(diffstrs, diffpans, tandotu))
    
    tanVels = diffstrs ./ diffpans

    pressCoeffs = [ pressure_coefficient((0, vt), margaret) for vt in tanVels ]
    
    liftCoeff = -sum(pressCoeffs .* (diffpans ./ 2) .* (cos.(panel_angle.(panels))))
    kuttaCoeff = -2 * strengths[end] / margaret

    return (liftCoeff, kuttaCoeff, pressCoeffs)
end

function parts(xs) 
    adj = collect(zip(drop(xs, 1), xs))
    return (first(adj), last(adj))
end 

adj2(xs) = collect(zip(drop(xs, 2), xs))

function midgrad(xs) 
    (firstpans, lastpans) = parts(xs)
    midpans = adj2(xs)
    return cat(dims = 1, firstpans, midpans, lastpans)
end

# Transforms (x, y) to the coordinate system with (x_s, y_s) as origin oriented at angle_s.
panel_coords(x, y, x_s, y_s, angle_s) = rotation(x - x_s, y - y_s, angle_s)
inverse_rotation(x, y, angle) = (x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle))
rotation(x, y, angle) = (x * cos(angle) + y * sin(angle), -x * sin(angle) + y * cos(angle))

pressure_coefficient(vels, mag) = 1 - (norm(vels))^2 / mag^2

# Performs velocity and potential calculations on a grid
function grid_data(objects::Array{<:Solution,1}, xs, pressure = true)
    vels = foldl((v1, v2)->[ u .+ v for (u, v) in zip(v1, v2)], [ velocity(object, xs) for object in objects ])
    pots = foldl((v1, v2)->v1 .+ v2, [ potential(object, xs) for object in objects ])
end

function grid_data(object::Solution, xs)
    vels = velocity(object, xs)
    pots = potential(object, xs)
    
    return vels, pots
end

# Performs velocity and potential computations for an object on a grid
velocity(object::Solution, xs) = [ velocity(object, x...) for x in xs ] 
potential(object::Solution, xs) = [ potential(object, x...) for x in xs ]

function cosine_panels(x::Array{<:Real}, y::Array{<:Real}, n = 40)
    """
    Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
    """
    d = (maximum(x) - minimum(x))
    x_center = (maximum(x) + minimum(x)) / 2.
    x_interp = cosineDist(x_center, d, n)

    itp = LinearInterpolation(x, y)
    y_interp = itp(x_interp)

    return (x_interp, y_interp)
end

function cosine_dist(x_center :: Real, diameter :: Real, n :: Integer = 40)
    """
    Provides the projections to the x-axis for a circle of given diameter and center.
    """
    return x_center .+ (diameter / 2.) .* cos.(range(-π, stop = 0, length = n))

function cosine_airfoil(x::Array{<:Real}, y::Array{<:Real}, n::Integer = 40)
    """
    Discretises a geometry consisting of x and y coordinates into panels by projecting the x-coordinate of a circle onto the geometry.
    """
    d = (maximum(x) - minimum(x))
    x_center = (maximum(x) + minimum(x)) / 2.
    x_circ = cosineDist(x_center, d, n)

    x_ends = copy(x_circ)
    y_ends = interpolate()

    return (x_ends, y_ends)
end

function NACA4(digits::Tuple, c, n, closed_te = false, split = false)
    
    # Airfoil characteristics
    # Camber
    m = digits[1] / 100
    # Position
    p = digits[2] / 10
    # Thickness-to-chord ratio
    t_by_c = (10 * digits[3] + digits[4]) / 100

    # Cosine spacing
    angles = range(0, stop = π, length = n + 1)
    xs = reverse([ c * (1 - 0.5 * (1 - cos(beta))) for beta in angles ])

    # Thickness distribution
    thickness = [ 5 * t_by_c * c * (0.2969 * sqrt(xc / c) - 0.126 * xc / c - 0.3516 * (xc / c)^2 + 0.2843 * (xc / c)^3 - (closed_te ? 0.1036 : 1015) * (xc / c)^4) for xc in xs ]
    
    if p == 0.0 || m == 0.0
        x_upper = xs
        y_upper = thickness
        x_lower = xs
        y_lower = -thickness
    else
        # Compute m
        mline = [ xc < p * c ? (m / p^2) * xc * (2 * p - xc / c) : (m * (c - xc) / (1 - p)^2) * (1 + xc / c - 2 * p) for xc in xs ]
        # Compute gradients
        gradients = [ xc < p * c ? atan((2 * m / p^2) * (p - xc / c)) : atan((2 * m / (1 - p)^2) * (p - xc / c)) for xc in xs ] 
        # Upper surface
        x_upper = xs .- thickness .* sin(gradients) 
        y_upper = mline .+ thickness .* cos(gradients)
        # Lower surface
        x_lower = xs .+ thickness .* sin(gradients) 
        y_lower = mline .- thickness .* cos(gradients)
    end

    if split
        upper = zip(x_upper, y_upper)
        lower = zip(x_lower, y_lower)
        return (reverse(upper), lower)
    else  
        (X, Y) = append!(reverse(x_upper), x_lower), append!(reverse(y_upper), y_lower)
        return (X, Y)
    end
end

end