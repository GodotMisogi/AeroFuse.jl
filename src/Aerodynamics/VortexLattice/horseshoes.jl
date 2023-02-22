## Horseshoe methods
#==========================================================================================#

function quarter_point(p1, p2) 
    μ = SVector(1/4, 0, 1/4)
    return @. (1 - μ) * p1 + μ * p2
end 

function three_quarter_point(p1, p2) 
    μ = SVector(3/4, 0, 3/4)
    return @. (1 - μ) * p1 + μ * p2
end 

control_point(p1, p2, p3, p4) = ( three_quarter_point(p1, p2) + three_quarter_point(p4, p3) ) / 2
bound_leg(p1, p2, p3, p4) = SVector(quarter_point(p1, p2), quarter_point(p4, p3))

"""
    bound_leg(panel :: Panel3D)

Compute the bound leg for a `Panel3D`, for horseshoes/vortex rings.
"""
bound_leg(panel :: Panel3D) = bound_leg(panel.p1, panel.p2, panel.p3, panel.p4)

"""
    control_point(panel :: Panel3D)

Compute the control point of a `Panel3D` for horseshoes/vortex rings, which is the 3-quarter point on each side in the ``x``-``z`` plane.
"""
control_point(panel :: Panel3D) = control_point(panel.p1, panel.p2, panel.p3, panel.p4)

# Velocity kernels
bound_leg_velocity(a, b, Γ)    = Γ/4π * (1/norm(a) + 1/norm(b)) * a × b / (norm(a) * norm(b) + dot(a, b))
trailing_leg_velocity(r, Γ, u) = Γ/4π * normalize(r) × normalize(u) / (norm(r) - dot(r, u))

trailing_legs_velocities(a, b, Γ, u) = trailing_leg_velocity(a, Γ, u) - trailing_leg_velocity(b, Γ, u)
total_horseshoe_velocity(a, b, Γ, u) = bound_leg_velocity(a, b, Γ) + trailing_legs_velocities(a, b, Γ, u)

# Finite-core velocity kernels
function bound_leg_velocity(a, b, Γ, ε) 
    na, nb, σ = norm(a), norm(b), dot(a, b)
    term_1 = (na^2 - σ) / √(na^2 + ε^2) + (nb^2 - σ) / √(nb^2 + ε^2)
    term_2 = a × b / (na^2 * nb^2 - σ^2 + ε^2 * (na^2 + nb^2 - 2 * na * nb))
    
    Γ/4π * term_1 * term_2
end

trailing_leg_velocity(r, Γ, u, ε) = Γ/4π * normalize(r) × u / (norm(r) - dot(r, u) + ε^2 / (norm(r) + dot(r, u)))
trailing_legs_velocities(a, b, Γ, u, ε) = trailing_leg_velocity(a, Γ, u, ε) - trailing_leg_velocity(b, Γ, u, ε)
total_horseshoe_velocity(a, b, Γ, u, ε) = bound_leg_velocity(a, b, Γ, ε) + trailing_legs_velocities(a, b, Γ, u, ε)

## Arrays of vortex lines
#==========================================================================================#

abstract type AbstractVortex end

"""
    Horseshoe(r1, r2, rc, normal, chord)

Define a horseshoe vortex with a start and endpoints ``r₁, r₂`` for the bound leg, a collocation point ``r``, a normal vector ``̂n``, and a finite core size.

The finite core setup is not implemented for now.
"""
struct Horseshoe{T <: Real} <: AbstractVortex
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
    rc :: SVector{3,T}
    normal :: SVector{3,T}
    core :: T
end

Base.length(::Horseshoe) = 1

r1(horseshoe :: Horseshoe) = horseshoe.r1
r2(horseshoe :: Horseshoe) = horseshoe.r2

function Horseshoe(r1, r2, rc, n, c)
    T = promote_type(eltype(r1), eltype(r2), eltype(rc), eltype(n), eltype(c))
    Horseshoe{T}(r1, r2, rc, n, c)
end

"""
    bound_leg(horseshoe :: Horseshoe)

Getter for bound leg field of a `Horseshoe`.
"""
bound_leg(horseshoe :: Horseshoe) = horseshoe.bound_leg
control_point(horseshoe :: Horseshoe)  = horseshoe.rc
normal_vector(horseshoe :: Horseshoe) = horseshoe.normal

r1(r, horseshoe :: Horseshoe) = r - horseshoe.r1
r2(r, horseshoe :: Horseshoe) = r - horseshoe.r2

"""
    Horseshoe(panel :: Panel3D, normal, drift = zeros(3))

Generate a `Horseshoe` corresponding to a `Panel3D`, an associated normal vector, and a "drift velocity".
"""
function Horseshoe(panel :: Panel3D, normal, drift = SVector(0., 0., 0.); core_size = 0.)
    r1, r2 = bound_leg(panel)
    rc = control_point(panel) + drift
    Horseshoe(r1, r2, rc, normal, core_size)
end

"""
    transform(horseshoe :: Horseshoe, T :: LinearMap)

Generate a new `Horseshoe` with the points and normal vectors transformed by the `LinearMap` ``T``.
"""
transform(horseshoe :: Horseshoe, T :: LinearMap) = setproperties(horseshoe,
    r1 = T(horseshoe.r1),
    r2 = T(horseshoe.r2),
    rc = T(horseshoe.rc),
    normal = T(horseshoe.normal),
)

transform(horseshoe :: Horseshoe; rotation = I(3), translation = zeros(3)) = transform(horseshoe, Translation(translation) ∘ LinearMap(rotation))

"""
    bound_leg_center(horseshoe :: Horseshoe)

Compute the midpoint of the bound leg of a `Horseshoe`.
"""
bound_leg_center(horseshoe :: Horseshoe) = (horseshoe.r1 + horseshoe.r2) / 2

"""
    bound_leg_vector(horseshoe :: Horseshoe)

Compute the direction vector of the bound leg of a `Horseshoe`.
"""
bound_leg_vector(horseshoe :: Horseshoe) = horseshoe.r2 - horseshoe.r1

"""
    velocity(r, horseshoe, Γ, u_hat = [1.,0.,0.])

Compute the induced velocity at a point ``r`` of a given `Horseshoe` with a bound leg of constant strength ``Γ`` and semi-infinite trailing legs pointing in a given direction ``û``, by default `û = x̂`.
"""
velocity(r, horseshoe :: Horseshoe, Γ :: Real, V_hat = SVector(1.,0.,0.)) = total_horseshoe_velocity(r1(r, horseshoe), r2(r, horseshoe), Γ, V_hat, horseshoe.core)

"""
    bound_velocity(r, horseshoe, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the bound leg with constant strength ``Γ`` of a given `Horseshoe`.
"""
bound_velocity(r, horseshoe :: Horseshoe, Γ) = bound_leg_velocity(r1(r, horseshoe), r2(r, horseshoe), Γ, horseshoe.core)

"""
    trailing_velocity(r, horseshoe, Γ, u_hat)

Compute the induced velocity at a point ``r`` from the semi-infinite trailing legs with constant strength ``Γ`` of a given `Horseshoe`.
"""
trailing_velocity(r, horseshoe :: Horseshoe, Γ, V) = trailing_legs_velocities(r1(r, horseshoe), r2(r, horseshoe), Γ, V, horseshoe.core)