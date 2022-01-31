## Reflections and projections
#==========================================================================================#

"""
    reflect_xz(vector)
    reflect_xz(line :: Line)
    reflect_xz(horseshoe :: Horseshoe)

Reflect the ``y``-coordinate of a given 3-dimensional vector (or ``Line`` or ``Horseshoe``) about the ``x``-``z`` plane.
"""
reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])
reflect_xz(line :: Line) = Line((reflect_xz ∘ p2)(line), (reflect_xz ∘ p1)(line))
reflect_xz(horseshoe :: Horseshoe) = (Horseshoe ∘ reflect_xz ∘ bound_leg)(horseshoe)

"""
    project_yz(vector)
    project_yz(line :: Line)

Project a given 3-dimensional vector or `Line` into the ``y``-``z`` plane.
"""
project_yz(vector) = SVector(0, vector[2], vector[3])
project_yz(line :: Line) = Line((project_yz ∘ p1)(line), (project_yz ∘ p2)(line))

"""
    flip_xz(vector)

Reflect the ``x``- and ``z``- coordinates of a given 3-dimensional vector about the ``y``-``z`` and ``x``-``y`` planes respectively for the representation in body axes.
"""
flip_xz(vector) = SVector(-vector[1], vector[2], -vector[3])

## Axis transformations
#==========================================================================================#

abstract type AbstractAxisSystem end

struct Geometry <: AbstractAxisSystem end
struct Body      <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind      <: AbstractAxisSystem end

Base.show(io :: IO, :: Geometry)  = print(io, "Geometry")
Base.show(io :: IO, :: Body)      = print(io, "Body")
Base.show(io :: IO, :: Stability) = print(io, "Stability")
Base.show(io :: IO, :: Wind)      = print(io, "Wind")


geometry_to_body_axes(coords) = flip_xz(coords)

"""
    geometry_to_stability_axes(coords, α)

Convert coordinates from geometry to stability axes with angle ``α``.
"""
geometry_to_stability_axes(coords, α :: T) where T <: Real = RotY{T}(α) * coords

"""
    geometry_to_stability_axes(coords, α)

Convert coordinates from stability to geometry axes with angle ``α``.
"""
stability_to_geometry_axes(coords, α :: T) where T <: Real = geometry_to_stability_axes(coords, -α)

# Possible cancellation errors causing convergence issues in optimization?
"""
    geometry_to_wind_axes(coords, α, β)
    geometry_to_wind_axes(line :: Line, α, β)

Convert coordinates from geometry axes to wind axes for given angles of attack ``α`` and sideslip ``\\beta.``
"""
geometry_to_wind_axes(coords, α, β) = let T = promote_type(eltype(α), eltype(β)); RotZY{T}(β, α) * coords end
# geometry_to_wind_axes(coords, α, β) = rotate_zy(β, α) * coords

# Define generated function?
geometry_to_wind_axes(line :: Line, α, β) = Line(geometry_to_wind_axes(r1(line), α, β), geometry_to_wind_axes(r2(line), α, β)) 

"""
    wind_to_geometry_axes(coords, α, β)

Convert coordinates from wind axes to geometry axes for given angles of attack ``α`` and sideslip \\beta.``
"""
wind_to_geometry_axes(coords, α, β) = let T = promote_type(eltype(α), eltype(β)); RotYZ{T}(-α, -β) * coords end

# function rotate_zy(θ₁, θ₂)
#     sinθ₁, cosθ₁ = sincos(θ₁)
#     sinθ₂, cosθ₂ = sincos(θ₂)
#     z            = zero(sinθ₁)

#     @SMatrix [ cosθ₁* cosθ₂  -sinθ₁  cosθ₁ * sinθ₂
#                sinθ₁* cosθ₂   cosθ₁  sinθ₁ * sinθ₂
#                  -sinθ₂         z        cosθ₂     ]
# end