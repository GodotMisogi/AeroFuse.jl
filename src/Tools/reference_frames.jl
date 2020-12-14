# Axis transformations
#==========================================================================================#

"""
Converts coordinates into stability axes.
"""
body_to_stability_axes(coords, α) = RotY{Float64}(α) * coords
                                            
"""
Converts coordinates from body axes to wind axes.
"""
body_to_wind_axes(coords, α, β) = RotZY{Float64}(β, α) * coords

"""
Reflects the y-coordinate of a given vector about the x-z plane.
"""
reflect_xz(vector :: SVector{3, <: Real}) = SVector(vector[1], -vector[2], vector[3])

"""
Reflects the x- and z- coordinates of a given vector about the y-z and x-y planes respectively.
"""
stab_flip(vector :: SVector{3, <: Real}) = SVector(-vector[1], vector[2], -vector[3])

"""
Transforms forces and moments from body to stability axes.
"""
body_to_stability_axes(force, moment, freestream :: Freestream) = body_to_stability_axes(force, freestream.α), body_to_stability_axes(stab_flip(moment), freestream.α), body_to_wind_axes(stab_flip(freestream.Ω), freestream.α)

"""
Transforms forces and moments from body to wind axes.
"""
body_to_wind_axes(force, moment, freestream :: Freestream) = body_to_wind_axes(force, freestream.α, freestream.β), body_to_wind_axes(stab_flip(moment), freestream.α, freestream.β), body_to_wind_axes(stab_flip(freestream.Ω), freestream.α, freestream.β)