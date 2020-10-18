module AircraftGeometry

include("MathTools.jl")
using LinearAlgebra
using .MathTools: structtolist

#----------------VECTOR SPACES?---------------#

abstract type Point end

*(scale :: Real, point :: Point) = typeof(point)(scale * structtolist(point)...)
+(p1 :: Point, p2 :: Point) = typeof(p1)(structtolist(p1) .+ structtolist(p2)...)

struct Point2D <: Point
    x :: Real; y :: Real;
end

struct Point3D <: Point
    x :: Real; y :: Real; z :: Real;
end

end