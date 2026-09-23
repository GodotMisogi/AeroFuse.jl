## Prototype for generic control deflection on panel mesh
#==========================================================================================#

abstract type AbstractControlSurface end

struct BoundingBox{T <: Real} <: AbstractControlSurface
    location :: SVector{3,T}
    length   :: T
    width    :: T
    function BoundingBox{T}(x_o, s_r, c_r) where T <: Real
        @assert all(0. < x < 1., x_o) "Origin of bounding box must lie within [0,1] ⊂ ℝ!"
        @assert 0. <= s_r <= 1. "The normalized span length must lie within [0,1] ⊂ ℝ!"
        @assert 0. <= c_r <= 1. "The normalized chord length must lie within [0,1] ⊂ ℝ!"
        
        new{T}(x_o, s_r, c_r)
    end
end

@Base.kwdef struct Flap{T} <: AbstractControlSurface
    angle = zero(T)
    hinge = convert(T, 0.75)
end

Flap(ang=0, hin=0.75) = Flap{promote_type(eltype(ang), eltype(hin))}(ang, hin)

@Base.kwdef struct Aileron{T} <: AbstractControlSurface
    angle = zero(T)
    hinge = convert(T, 0.75)
end

Aileron(ang=0, hin=0.75)  = Aileron{promote_type(eltype(ang), eltype(hin))}(ang, hin)

const Elevator = Flap
const Rudder = Flap