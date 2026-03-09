## Reflections and projections
#==========================================================================================#

# Reflect the ``y``-coordinate of a given 3-dimensional vector about the ``x``-``z`` plane.
reflect_xz(vector) = SVector(vector[1], -vector[2], vector[3])

# Project a given 3-dimensional vector or into the ``y``-``z`` plane.
project_yz(vector) = SVector(0, vector[2], vector[3])

# Reflect the ``x``- and ``z``- coordinates of a given 3-dimensional vector about the ``y``-``z`` and ``x``-``y`` planes respectively for the representation in body axes.
flip_xz(vector) = SVector(-vector[1], vector[2], -vector[3])

## Axis transformations
#==========================================================================================#

# Convert coordinates from geometry to body axes
geometry_to_body_axes(xyz) = flip_xz(xyz)

# Convert coordinates from geometry to stability axes with angle ``α``.
geometry_to_stability_axes(xyz, α :: T) where T <: Real = RotY{T}(α) * xyz

# Convert coordinates from stability to geometry axes with angle ``α``.
stability_to_geometry_axes(xyz, α :: T) where T <: Real = geometry_to_stability_axes(xyz, -α)

# Convert coordinates from geometry axes to wind axes for given angles of attack ``α`` and sideslip ``β``.
geometry_to_wind_axes(xyz, α, β) = let T = promote_type(eltype(xyz), eltype(α), eltype(β)); RotZY{T}(β, α) * xyz end

## Check order
function wind_to_geometry_axes(xyz, α, β) 
    T = promote_type(eltype(α), eltype(β))
    return RotYZ{T}(α, β) * xyz
end

wind_to_body_axes(xyz, α, β) = flip_xz(wind_to_geometry_axes(xyz, α, β))

## Generic axis transformation helpers
#==========================================================================================#

# Transform a vector (force/velocity) from geometry axes to the target axis system.
_to_axes(v, ::Geometry, α, β) = v
_to_axes(v, ::Body, α, β) = geometry_to_body_axes(v)
_to_axes(v, ::Stability, α, β) = geometry_to_stability_axes(v, α)
_to_axes(v, ::Wind, α, β) = geometry_to_wind_axes(v, α, β)

# Transform a moment vector from geometry axes to the target axis system, applying the sign flip convention for stability and wind axes.
_moment_to_axes(m, ::Geometry, α, β) = m
_moment_to_axes(m, ::Body, α, β) = geometry_to_body_axes(m)
_moment_to_axes(m, ::Stability, α, β) = geometry_to_stability_axes(flip_xz(m), α)
_moment_to_axes(m, ::Wind, α, β) = geometry_to_wind_axes(flip_xz(m), α, β)