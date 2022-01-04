## Horseshoe methods
#==========================================================================================#

quarter_point(p1, p2) = weighted_vector(p1, p2, SVector(1/4, 0, 1/4))

three_quarter_point(p1, p2) = weighted_vector(p1, p2, SVector(3/4, 0, 3/4))

collocation_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) + three_quarter_point(p4, p3) ) / 2
bound_leg(p1, p2, p3, p4) = SVector(quarter_point(p1, p2), quarter_point(p4, p3))

"""
    bound_leg(panel :: Panel3D)

Compute the bound leg for a `Panel3D`, for horseshoes/vortex rings.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
    collocation_point(panel :: Panel3D)

Compute the collocation point of a `Panel3D` for horseshoes/vortex rings, which is the 3-quarter point on each side in the ``x``-``z`` plane.
"""
horseshoe_point(panel :: Panel3D) = collocation_point(panel.p1, panel.p2, panel.p3, panel.p4)

abstract type AbstractLine end

"""
    Line(r1 :: SVector{3,<: Real}, r2 :: SVector{3,<: Real})

Define a line segment based on two 3-dimensional vectors.
"""
struct Line{T <: Real} <: AbstractLine
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
end

Line(r1, r2) = let T = promote_type(eltype(r1), eltype(r2)); Line{T}(r1, r2) end
Line((r1, r2)) = Line(r1, r2)

r1(line :: Line) = line.r1
r2(line :: Line) = line.r2
vector(line :: Line) = r2(line) - r1(line)
center(line :: Line) = (r1(line) + r2(line)) / 2

points(lines) = [ r1.(lines); [(r2 ∘ last)(lines)] ]

transform(line :: Line, rotation, translation) = let trans = Translation(translation) ∘ LinearMap(rotation); Line((trans ∘ r1)(line), (trans ∘ r2)(line)) end

r1(r, line :: Line) = r - r1(line)
r2(r, line :: Line) = r - r2(line)

bound_leg_velocity(a, b, Γ)    = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))
trailing_leg_velocity(r, Γ, u) = Γ/4π * normalize(r) × normalize(u) / (norm(r) - dot(r, u))

trailing_legs_velocities(a, b, Γ, u) = trailing_leg_velocity(a, Γ, u) - trailing_leg_velocity(b, Γ, u)
total_horseshoe_velocity(a, b, Γ, u) = bound_leg_velocity(a, b, Γ) + trailing_legs_velocities(a, b, Γ, u)

horseshoe_velocity(r, line :: Line, Γ, direction) = total_horseshoe_velocity(r1(r, line), r2(r, line), Γ, direction)

bound_velocity(r, line :: Line, Γ) = bound_leg_velocity(r1(r, line), r2(r, line), Γ)

## Arrays of vortex lines
#==========================================================================================#

abstract type AbstractVortexArray end

"""
    Horseshoe(bound_leg, collocation_point, normal, chord)

Define a horseshoe vortex based on a bound leg, a collocation point, a normal vector, and a "chord length".
"""
struct Horseshoe{T <: Real} <: AbstractVortexArray
    bound_leg         :: Line{T}
    collocation_point :: SVector{3,T}
    normal            :: SVector{3,T}
    chord             :: T
end

Base.length(::Horseshoe) = 1

"""
    bound_leg(horseshoe :: Horseshoe)

Getter for bound leg field of a `Horseshoe`.
"""
bound_leg(horseshoe :: Horseshoe)        = horseshoe.bound_leg
horseshoe_point(horseshoe :: Horseshoe)  = horseshoe.collocation_point
horseshoe_normal(horseshoe :: Horseshoe) = horseshoe.normal

r1(r, horseshoe :: Horseshoe) = r1(r, bound_leg(horseshoe))
r2(r, horseshoe :: Horseshoe) = r2(r, bound_leg(horseshoe))

"""
    Horseshoe(panel :: Panel3D, normal, drift = zeros(3))

Generate a `Horseshoe` corresponding to a `Panel3D`, an associated normal vector, and a "drift velocity".
"""
Horseshoe(panel :: Panel3D, normal, drift = SVector(0., 0., 0.)) =
    Horseshoe(Line(bound_leg(panel)), horseshoe_point(panel) + drift, normal, (norm ∘ average_chord)(panel))

"""
Compute the midpoint of the bound leg of a `Horseshoe`.
"""
bound_leg_center(horseshoe :: Horseshoe) = (center ∘ bound_leg)(horseshoe)

"""
Compute the direction vector of the bound leg of a `Horseshoe`.
"""
bound_leg_vector(horseshoe :: Horseshoe) = (vector ∘ bound_leg)(horseshoe)

"""
    velocity(r, horseshoe, Γ, V_hat)

Compute the induced velocities at a point ``r`` of a given Horseshoe with constant strength ``Γ`` and trailing legs pointing in a given direction ``û``.
"""
function velocity(r, horseshoe :: Horseshoe, Γ :: Real, V_hat, finite_core = false) 
    if finite_core
        width = (norm ∘ bound_leg_vector)(horseshoe)
        ε = 0.1 * max(horseshoe.chord, width) # Wrong core size? Consider options...
        horseshoe_velocity(r, bound_leg(horseshoe), Γ, V_hat, ε)
    else
        horseshoe_velocity(r, bound_leg(horseshoe), Γ, V_hat)
    end
end

"""
    trailing_velocity(r, horseshoe, Γ, V_hat)

Compute the induced velocities at a point ``r`` of a given Horseshoe with constant strength ``Γ`` and trailing legs pointing in a given direction ``û``, excluding the bound leg contribution.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r1(r, horseshoe), r2(r, horseshoe), Γ, V)

trailing_velocity(r, horseshoe :: Horseshoe, Γ, V, ε) = trailing_legs_velocities(r1(r, horseshoe), r2(r, horseshoe), Γ, V, ε)