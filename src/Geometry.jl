module Geometry

#----------------VECTOR SPACES?---------------#

abstract type Point end

struct Point2D <: FieldVector{2, Float64}
    x :: Real; y :: Real;
end

struct Point3D <: FieldVector{3, Float64}
    x :: Real; y :: Real; z :: Real;
end

yflip!(xs) = xs[:,2] .= -xs[:,2]

end