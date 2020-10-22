module Geometry

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

end