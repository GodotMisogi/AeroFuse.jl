## Standard stabilizers (non-exhaustive)
#==========================================================================================#

HorizontalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_twist = 0., tip_twist = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil, angle = 0.) = WingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, root_twist = root_twist, tip_twist = tip_twist, root_foil = naca4(0,0,1,2), tip_foil = tip_foil, angle = angle, axis = [0., 1., 0.])

VerticalTail(; root_chord = 1.0, taper = 1.0, span = 1.0, sweep = 0., w_sweep = 0., root_foil = naca4(0,0,1,2), tip_foil = root_foil) = HalfWingSection(root_chord = root_chord, taper = taper, span = span, sweep = sweep, w_sweep = w_sweep, tip_foil = tip_foil, angle = 90., axis = [1., 0., 0.])

## Prototype for generic control deflection on panel mesh
#==========================================================================================#

abstract type AbstractControlSurface end

struct BoundingBox{T <: Real} <: AbstractControlSurface
    location    :: SVector{3,T}
    span_ratio  :: T
    chord_ratio :: T
    function BoundingBox{T}(x_o, s_r, c_r) where T <: Real
        @assert all(0. < x < 1., x_o) "Origin of bounding box must lie within [0,1] ⊂ ℝ!"
        @assert 0. <= s_r <= 1. "The normalized span length must lie within [0,1] ⊂ ℝ!"
        @assert 0. <= c_r <= 1. "The normalized chord length must lie within [0,1] ⊂ ℝ!"
        
        new{T}(x_o, s_r, c_r)
    end
end

function Aileron(bb :: BoundingBox, δ)
    
end