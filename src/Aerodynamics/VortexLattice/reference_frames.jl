## Axis transformations
#==========================================================================================#

abstract type AbstractAxisSystem end

struct Global    <: AbstractAxisSystem end
struct Body      <: AbstractAxisSystem end
struct Stability <: AbstractAxisSystem end
struct Wind      <: AbstractAxisSystem end


function rotate_zy(theta1, theta2)
    sinθ₁, cosθ₁ = sincos(theta1)
    sinθ₂, cosθ₂ = sincos(theta2)
    z = zero(sinθ₁)

    @SMatrix [ cosθ₁* cosθ₂  -sinθ₁  cosθ₁ * sinθ₂
               sinθ₁* cosθ₂   cosθ₁  sinθ₁ * sinθ₂
                 -sinθ₂         z        cosθ₂     ]
end


"""
    body_to_stability_axes(coords, α)

Convert coordinates into stability axes with angle ``\\alpha``.
"""
body_to_stability_axes(coords, α :: T) where T <: Real = RotY{T}(α) * coords

"""
    body_to_stability_axes(coords, α)

Convert coordinates into stability axes with angle ``\\alpha``.
"""
stability_to_body_axes(coords, α :: T) where T <: Real = RotY{T}(-α) * coords
                                            

# Possible cancellation errors causing convergence issues in optimization?
# body_to_wind_axes(coords, α, β) = let T = promote_type(eltype(α), eltype(β)); RotZY{T}(β, α) * coords end

"""
    body_to_wind_axes(coords, α, β)

Convert coordinates from body axes to wind axes with angles ``\\alpha,~ \\beta``.
"""
body_to_wind_axes(coords, α, β) = rotate_zy(β, α) * coords

"""
    body_to_wind_axes(coords, α, β)

Convert coordinates from wind axes to body axes with angles ``\\alpha,~ \\beta``.
"""
wind_to_body_axes(coords, α :: T, β :: T) where T <: Real = RotZY{T}(-α, -β) * coords

## Line methods
#==========================================================================================#

"""
    body_to_wind_axes(line, freestream)

Transform a Line from body to wind axes in a given `Freestream`.
"""
body_to_wind_axes(line :: Line, α, β) = Line(body_to_wind_axes(r1(line), α, β), body_to_wind_axes(r2(line), α, β)) 


## Reflections and projections
#==========================================================================================#

"""
    reflect_xz(vector)

Reflect the y-coordinate of a given 3-dimensional vector about the x-z plane.
"""
reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])

"""
    project_yz(vector)

Project a given 3-dimensional vector into the y-z plane.
"""
project_yz(vector) = SVector(0, vector[2], vector[3])

"""
    reflect_yz(line)

Reflect a Line onto the ``x``-``z`` plane of its reference coordinate system.
"""
reflect_xz(line :: Line) = Line((reflect_xz ∘ p2)(line), (reflect_xz ∘ p1)(line))

"""
    reflect_yz(horseshoe)

Reflect a Horseshoe onto the ``x``-``z`` plane of its reference coordinate system.
"""
reflect_xz(horseshoe :: Horseshoe) = (Horseshoe ∘ reflect_xz ∘ bound_leg)(horseshoe)

"""
    reflect_yz(panel :: Panel3D)

Reflect a Panel3D onto the ``x``-``z`` plane of its reference coordinate system.
"""
reflect_xz(panel :: Panel3D) = Panel3D((reflect_xz ∘ p1)(panel), (reflect_xz ∘ p2)(panel), (reflect_xz ∘ p3)(panel), (reflect_xz ∘ p4)(panel))

"""
    project_yz(line)

Project a Line onto the ``y``-``z`` plane of its reference coordinate system.
"""
project_yz(line :: Line) = Line((project_yz ∘ p1)(line), (project_yz ∘ p2)(line))

"""
    stability_flip(vector)

Reflect the x- and z- coordinates of a given 3-dimensional vector about the y-z and x-y planes respectively.
"""
stability_flip(vector) = SVector(-vector[1], vector[2], -vector[3])