# Axis transformations
#==========================================================================================#

"""
    body_to_stability_axes(coords, α)

Converts coordinates into stability axes with angle ``\\alpha``.
"""
body_to_stability_axes(coords, α :: Real) = RotY{Float64}(α) * coords
                                            
"""
    body_to_wind_axes(coords, α, β)

Converts coordinates from body axes to wind axes with angles ``\\alpha,~ \\beta``.
"""
body_to_wind_axes(coords, α :: Real, β :: Real) = RotZY{Float64}(β, α) * coords

"""
    body_to_wind_axes(coords, α, β)

Converts coordinates from body axes to wind axes with angles ``\\alpha,~ \\beta``.
"""
body_to_wind_axes(vector :: SVector{3, <: Real}, freestream :: Freestream) = body_to_wind_axes(vector, freestream.α, freestream.β)

"""
    reflect_xz(vector)

Reflects the y-coordinate of a given SVector about the x-z plane.
"""
reflect_xz(vector :: SVector{3, <: Real}) = SVector(vector[1], -vector[2], vector[3])

"""
    project_yz(line)

Projects a Line onto the y-z plane of its coordinate system.
"""
project_yz(line :: Line) = Line(SVector(0, line.r1[2], line.r1[3]), SVector(0, line.r2[2], line.r2[3]))

"""
    stability_flip(vector)

Reflects the x- and z- coordinates of a given SVector about the y-z and x-y planes respectively.
"""
stability_flip(vector :: SVector{3, <: Real}) = SVector(-vector[1], vector[2], -vector[3])

"""
    body_to_wind_axes(line, freestream)

Transforms a Line from body to wind axes in a given Freestream.
"""
body_to_wind_axes(line :: Line, freestream :: Freestream) = Line(body_to_wind_axes(line.r1, freestream.α, freestream.β), body_to_wind_axes(line.r2, freestream.α, freestream.β)) 

"""
    body_to_stability_axes(force, moment, freestream)

Transforms forces and moments from body to stability axes in a given Freestream.
"""
body_to_stability_axes(force, moment, freestream :: Freestream) = body_to_stability_axes(force, freestream.α), body_to_stability_axes(stability_flip(moment), freestream.α), body_to_stability_axes(stability_flip(freestream.Ω), freestream.α)

"""
    body_to_wind_axes(force, moment, freestream)

Transforms forces and moments from body to wind axes in a given Freestream.
"""
body_to_wind_axes(force, moment, freestream :: Freestream) = body_to_wind_axes(force, freestream.α, freestream.β), body_to_wind_axes(stability_flip(moment), freestream.α, freestream.β), body_to_wind_axes(stability_flip(freestream.Ω), freestream.α, freestream.β)