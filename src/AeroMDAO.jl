module AeroMDAO

abstract type Panel3D end

struct DoubletPanel3D <: Panel3D
    p1 :: Tuple{Float64, Float64, Float64}
    p2 :: Tuple{Float64, Float64, Float64}
    p3 :: Tuple{Float64, Float64, Float64}
    p4 :: Tuple{Float64, Float64, Float64}
end

end # module
