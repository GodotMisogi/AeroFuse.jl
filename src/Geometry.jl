module Geometry

export Panel, Panel2D, Panel3D, split_panels, panel_angle, panel_length, panel_normal, panel_tangent, panel_dist, collocation_point, make_panels

include("MathTools.jl")
include("LaplaceSolutions.jl")
using LinearAlgebra
using StaticArrays
using .MathTools: span, structtolist

#----------------VECTOR SPACES?---------------#

abstract type Point end

struct Point2D <: FieldVector{2, Float64}
    x :: Real; y :: Real;
end

struct Point3D <: FieldVector{3, Float64}
    x :: Real; y :: Real; z :: Real;
end

abstract type Panel <: Laplace end

# Methods on panels in N dimensions
collocation_point(panel :: Panel) = let panel = structtolist(panel); sum(panel) / len(panel) end
panel_dist(panel_1 :: Panel, panel_2 :: Panel) = norm(collocation_point(panel_2) .- collocation_point(panel_1))
split_panels(panels :: Array{<: Panel}) = collect.(span(panel -> panel_location(panel) == "upper", panels))

struct Panel2D <: Panel
    p1 :: Point2D
    p2 :: Point2D
end

make_panels(coords :: Array{<: Real, 2}) = [ Panel2D(Point2D(xs, ys), Point2D(xe, ye)) for (xs, ys, xe, ye) ∈ eachrow([ coords[2:end,:] coords[1:end-1,:] ]) ][end:-1:1]

panel_length(panel :: Panel2D) = norm(panel.p2 .- panel.p1)
panel_angle(panel :: Panel2D) = let (xs, ys) = panel.p1, (xe, ye) = panel.p2; atan(ye - ys, xe - xs) end
panel_tangent(panel :: Panel2D) = rotation(1, 0, -1 * panel_angle(panel))
panel_normal(panel :: Panel2D) = inverse_rotation(0, 1, panel_angle(panel))
panel_location(panel :: Panel2D) = let angle = panel_angle(panel); (π/2 <= angle <= π) || (-π <= angle <= -π/2) ? "lower" : "upper" end

struct Panel3D <: Panel
    p1 :: Point3D
    p2 :: Point3D
    p3 :: Point3D
    p4 :: Point3D
end

panel_normal(panel :: Panel3D) = cross(panel.p2 .- panel.p1, panel.p3 .- panel.p2)

struct Line
    r1 :: Point
    r2 :: Point
end

end