struct WingCS{T <: Real} <: Aircraft
    wing :: Union{Wing{T}, HalfWing{T}}
    position :: SVector{3,T}
    orientation :: Matrix{T}
end

mean_aerodynamic_center(comp :: WingCS) = (Translation(comp.position) ∘ LinearMap(comp.orientation)) (mean_aerodynamic_center(comp.wing))