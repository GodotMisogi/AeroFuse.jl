struct Coordinates{T <: Real} <: Aircraft
    position :: SVector{3,T}
    axis     :: SVector{3,T}
    angle    :: T
end

wing(comp :: WingCS) = comp.wing
position(comp :: WingCS) = comp.position
axis(comp :: WingCS) = comp.axis
angle(comp :: WingCS) = comp.angle
orientation(comp :: WingCS{T}) where T <: Real = AngleAxis{T}(comp.angle, comp.axis...)

affine_transformation(comp :: WingCS) = Translation(position(comp)) ∘ LinearMap(orientation(comp))

mean_aerodynamic_center(trans, wing :: Union{HalfWing, Wing}) = affine_transformation(trans)(mean_aerodynamic_center(wing))