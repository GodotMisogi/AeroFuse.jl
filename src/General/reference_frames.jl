# Axis transformations
#==========================================================================================#

"""
Converts coordinates into stability axes.
"""
body_to_stability_axes(coords, uni :: Freestream) =  
                                            RotY{Float64}(uni.α) * coords
                                            # [cos(uni.α) 0 -sin(uni.α); 
                                            #       0     1     0     ;
                                            #  sin(uni.α) 0 sin(uni.α)] * coords
                                            

"""
Converts coordinates from body axes to wind axes.
"""
body_to_wind_axes(coords, uni :: Freestream) = RotZY{Float64}(uni.β, uni.α) * coords

"""
Reflects the y-coordinate of a given vector about the x-z plane.
"""
reflect_xz(vector :: SVector{3, Float64}) = SVector(vector[1], -vector[2], vector[3])

"""
Reflects the x- and z- coordinates of a given vector about the y-z and x-y planes respectively.
"""
stab_flip(vector :: SVector{3, Float64}) = SVector(-vector[1], vector[2], -vector[3])

"""
Transforms forces and moments from body to stability axes.
"""
body_to_stability_axes(force, moment, Ω, freestream :: Freestream) = body_to_stability_axes(force, freestream), body_to_stability_axes(stab_flip(moment), freestream), body_to_stability_axes(stab_flip(Ω), freestream)

"""
Transforms forces and moments from body to wind axes.
"""
body_to_wind_axes(force, moment, Ω, freestream :: Freestream) = body_to_wind_axes(force, freestream), body_to_wind_axes(stab_flip(moment), freestream), body_to_wind_axes(stab_flip(Ω), freestream)

"""
Transforms forces and moments into wind axes.
"""
body_to_wind_axes(force, moment, freestream :: Freestream) = body_to_wind_axes(force, freestream), body_to_wind_axes([ -moment[1], moment[2], -moment[3] ], freestream)
