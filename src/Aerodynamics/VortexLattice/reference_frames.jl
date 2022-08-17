## Reflections and projections
#==========================================================================================#

"""
    reflect_xz(vector)
    reflect_xz(line :: Line)
    reflect_xz(horseshoe :: Horseshoe)

Reflect the ``y``-coordinate of a given 3-dimensional vector (or ``Line`` or ``Horseshoe``) about the ``x``-``z`` plane.
"""
reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])
reflect_xz(horseshoe :: Horseshoe) = (Horseshoe ∘ reflect_xz ∘ bound_leg)(horseshoe)

"""
    project_yz(vector)
    project_yz(line :: Line)

Project a given 3-dimensional vector or `Line` into the ``y``-``z`` plane.
"""
project_yz(vector) = SVector(0, vector[2], vector[3])

"""
    flip_xz(vector)

Reflect the ``x``- and ``z``- coordinates of a given 3-dimensional vector about the ``y``-``z`` and ``x``-``y`` planes respectively for the representation in body axes.
"""
flip_xz(vector) = SVector(-vector[1], vector[2], -vector[3])

## Axis transformations
#==========================================================================================#

geometry_to_body_axes(coords) = -coords

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

"""
    geometry_to_wind_axes(coords, α, β)
    geometry_to_wind_axes(horseshoe :: Horseshoe, α, β)

Convert coordinates from geometry axes to wind axes for given angles of attack ``α`` and sideslip ``\\beta.``
"""
geometry_to_wind_axes(coords, α, β) = let T = promote_type(eltype(α), eltype(β)); RotZY{T}(β, α) * coords end

function geometry_to_wind_axes(horseshoe :: Horseshoe, α, β) 
    T = promote_type(eltype(α), eltype(β))
    return transform(horseshoe, LinearMap(RotZY{T}(β, α)))
end

geometry_to_wind_axes(coords, fs :: Freestream) = geometry_to_wind_axes(coords, fs.alpha, fs.beta)
geometry_to_wind_axes(horseshoe :: Horseshoe, fs :: Freestream) = geometry_to_wind_axes(horseshoe, fs.alpha, fs.beta)

"""
    wind_to_geometry_axes(coords, α, β)
    wind_to_geometry_axes(horseshoe :: Horseshoe, α, β) 

Convert coordinates from wind axes to geometry axes for given angles of attack ``α`` and sideslip \\beta.``
"""
## Check order
function wind_to_geometry_axes(coords, α, β) 
    T = promote_type(eltype(α), eltype(β))
    return RotYZ{T}(-α, -β) * coords
end

function wind_to_geometry_axes(horseshoe :: Horseshoe, α, β) 
    T = promote_type(eltype(α), eltype(β))
    return transform(horseshoe, LinearMap(RotYZ{T}(-α, -β)))
end


# function rotate_zy(θ₁, θ₂)
#     sinθ₁, cosθ₁ = sincos(θ₁)
#     sinθ₂, cosθ₂ = sincos(θ₂)
#     z            = zero(sinθ₁)

#     @SMatrix [ cosθ₁* cosθ₂  -sinθ₁  cosθ₁ * sinθ₂
#                sinθ₁* cosθ₂   cosθ₁  sinθ₁ * sinθ₂
#                  -sinθ₂         z        cosθ₂     ]
# end