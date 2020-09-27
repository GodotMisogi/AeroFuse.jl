module AeroMDAO

using LinearAlgebra
using Base.Iterators
using Interpolations
using StaticArrays

abstract type Panel end

abstract type Aircraft end

struct Panel3D <: Panel
    loc :: Array{Tuple{Float64, Float64, Float64}, 4}
end

struct Panel2D <: Panel
    loc :: Array{Tuple{Float64, Float64}, 2}
end

# Methods on panels
midpoint(panel :: Panel) = mean(panel.loc, 1)

struct WingSection <: Aircraft
    """
    Definition for a wing cross section ("X-section").
    """
    
    loc :: Array{Tuple{Float64, Float64, Float64}} # Coordinates of the section's leading edge with respect to the reference frame of the wing.
    chord :: Float64  # Chord length of the section
    twist :: Float64  # Twist of the section about the leading edge
    airfoil :: Array{Tuple{Float64, Float64}} # The airfoil profile as an array of coordinates
end

end
