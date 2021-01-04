"""
Helper function to compute the bound leg of Panel3D for horseshoes/vortex rings, which is the quarter point on each side in the x-z plane.
"""
bound_leg(p1, p2, p3, p4) = [ quarter_point(p1, p2), quarter_point(p4, p3) ]

"""
Helper function to compute the collocation point of Panel3D for horseshoes/vortex rings, which is the 3-quarter point on each side in the x-z plane.
"""
collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) + three_quarter_point(p4, p3) ) / 2

"""
Computes the bound leg for a Panel3D.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
Computes the collocation point of a Panel3D.
"""
collocation_point(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

"""
A composite type consisting of two vectors to define a line.
"""
struct Line
    r1 :: SVector{3, Real}
    r2 :: SVector{3, Real}
end

point1(line :: Line) = line.r1
point2(line :: Line) = line.r2
vector(line :: Line) = point2(line) - point1(line)
center(line :: Line) = (point1(line) + point2(line)) / 2

points(lines :: AbstractVector{Line}) = [ point1.(lines)..., (point2 ∘ last)(lines) ]

transform(line :: Line, rotation, translation) = let trans = Translation(translation) ∘ LinearMap(rotation); Line((trans ∘ point1)(line), (trans ∘ point2)(line)) end

"""
Helper function to compute the velocity induced by a bound vortex leg.
"""
bound_leg_velocity(a, b, Γ) = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))

"""
Helper function to compute the velocity induced by trailing vortex legs.
"""
trailing_legs_velocities(a, b, Γ, u) = Γ/4π * (a × u / (norm(a) - dot(a, u)) / norm(a) - b × u / (norm(b) - dot(b, u)) / norm(b))

"""
Placeholder.
"""
total_horseshoe_velocity(a, b, Γ, u) = bound_leg_velocity(a, b, Γ) + trailing_legs_velocities(a, b, Γ, u)

"""
Computes the velocity induced at a point `r` by a vortex Line with constant strength Γ.
"""
horseshoe_velocity(r, line :: Line, Γ, direction) = @timeit "Sum Legs" total_horseshoe_velocity(r - point1(line), r - point2(line), Γ, direction)
## Arrays of vortex lines
#==========================================================================================#

"""
Placeholder. Vortex rings and horseshoes basically have the same methods, and are arrays of vortex lines.
"""
abstract type AbstractVortexArray end

"""
A horseshoe type consisting of a bound leg of type Line represening a vortex line.
"""
struct Horseshoe <: AbstractVortexArray
    bound_leg :: Line
end

"""
    bound_leg(horseshoe)

Getter for bound leg field of a Horseshoe.
"""
bound_leg(horseshoe :: Horseshoe) = horseshoe.bound_leg

"""
Returns a Horseshoe on a Panel3D.
"""
horseshoe_lines(panel :: Panel3D) = (Horseshoe ∘ Line)(bound_leg(panel)...)

"""
Computes the midpoint of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_center(horseshoe :: Horseshoe) = (center ∘ bound_leg)(horseshoe)

"""
Computes the direction vector of the bound leg of a Horseshoe or Vortex Ring.
"""
bound_leg_vector(horseshoe :: Horseshoe) = (vector ∘ bound_leg)(horseshoe)

"""
    velocity(r, horseshoe, Γ, V_hat)

Computes the induced velocities at a point ``r`` of a given Horseshoe with constant strength ``Γ`` and trailing legs pointing in a given direction ``\\hat V``.
"""
velocity(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, Γ :: Real, V_hat :: SVector{3, <: Real}) = horseshoe_velocity(r, bound_leg(horseshoe), Γ, V_hat)

"""
Placeholder.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = let line = bound_leg(horseshoe); @timeit "Trailing Velocity" trailing_legs_velocities(r - point1(line), r - point2(line), Γ, V) end