using ForwardDiff, ReverseDiff, Zygote
using LinearAlgebra, StaticArrays, Rotations

## Foil type
#===================================================================#

abstract type AbstractFoil end

struct Foil{T <: Real} <: AbstractFoil
    x :: Vector{T}
    y :: Vector{T}
end

Foil(x :: AbstractVector{T}, y :: AbstractVector{T}) where T <: Real = Foil{T}(x, y)
Foil(coords) = Foil(coords[:,1], coords[:,2])

coordinates(foil :: Foil) = [ foil.x foil.y ]

arc_length(foil :: Foil) = let c = coordinates(foil); norm(c[2:end,:] .- c[1:end-1,:]) end

# Coordinates
x_coords = [ 1.0  0.0
             0.5  0.5
             0.0  0.0
             0.5 -0.5
             1.0  0.0 ]

## ForwardDiff (PASSES)
ForwardDiff.gradient(arc_length ∘ Foil, x_coords)

## ReverseDiff (PASSES)
ReverseDiff.gradient(arc_length ∘ Foil, x_coords)

## Zygote (PASSES)
Zygote.gradient(arc_length ∘ Foil, x_coords)
Zygote.gradient(arc_length ∘ Foil, x_coords[:,1], x_coords[:,2])

## Bing type
#===================================================================#

struct Bing{T <: Real, N <: AbstractFoil}
    foils  :: Vector{N}
    chords :: Vector{T}
end

Bing(foils :: AbstractVector{N}, chords :: AbstractVector{T}) where {T <: Real, N <: AbstractFoil} = Bing{T,N}(foils, chords)

arc_length(fw :: Bing) = sum(arc_length, fw.foils)

foiler_wing(x1, x2) = arc_length(Bing(fill(Foil(x1), length(x2)), x2))

diff_foiler_wing(x) = foiler_wing(reshape(x[1:end-length(cs)], size(x_coords)), x[end-length(cs):end])

##
cs  = [1.0, 0.6, 0.3]

foilwing = foiler_wing(x_coords, cs)

##
diff_foiler_wing([ x_coords[:]; cs ])

## ForwardDiff (PASSES)
ForwardDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs ])  
ForwardDiff.gradient(x -> arc_length(Bing(fill(Foil(x), length(cs)), cs)), x_coords)
ForwardDiff.gradient(x -> arc_length(Bing(fill(Foil(x_coords), length(x)), x)), cs)

## ReverseDiff (PASSES)
ReverseDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs ])
ReverseDiff.gradient(x -> arc_length(Bing(fill(Foil(x), length(cs)), cs)), x_coords)
ReverseDiff.gradient(x -> arc_length(Bing(fill(Foil(x_coords), length(x)), x)), cs)

## Zygote (FAILS)
Zygote.gradient(diff_foiler_wing, [ x_coords[:]; cs ])

## Line type
#===================================================================#

abstract type AbstractLine end

struct Line{T <: Real} <: AbstractLine
    r1 :: SVector{3,T}
    r2 :: SVector{3,T}
end

function rotate_zy(theta1, theta2)
    sinθ₁, cosθ₁ = sincos(theta1)
    sinθ₂, cosθ₂ = sincos(theta2)
    z = zero(sinθ₁)

    # transposed representation
    [ cosθ₁*cosθ₂  sinθ₁*cosθ₂ -sinθ₂
     -sinθ₁        cosθ₁        z
      cosθ₁*sinθ₂  sinθ₁*sinθ₂  cosθ₂ ]
end

Line(r1 :: AbstractVector{M}, r2 :: AbstractVector{N}) where {M <: Real, N <: Real} = Line{M,N}(r1, r2)
Line((r1, r2)) = Line(r1, r2)

body_to_wind_axes(coords, α :: T, β :: T) where T <: Real = rotate_zy(β, α) * coords
body_to_wind_axes(line :: AbstractLine, α :: T, β :: T) where {T <: Real} = Line(body_to_wind_axes(line.r1, α, β), body_to_wind_axes(line.r2, α, β))

##
α  = deg2rad(1)
β  = deg2rad(2)
r1 = [1,2.,3]
r2 = [1,0.,1]
x  = [ r1 r2 ]
line = Line(x[:,1], x[:,2])

## ForwardDiff
ForwardDiff.gradient(x -> sum(body_to_wind_axes(x - r2, α, β)), r1)
ForwardDiff.gradient(x -> let l = Line(x[:,1], x[:,2]); sum([ l.r1; l.r2 ]) end, x)
ForwardDiff.gradient(x -> let l = Line(x[:,1], x[:,2]); sum(body_to_wind_axes(sum([ l.r1 l.r2 ], dims = 2), α, β)) end, x)
ForwardDiff.gradient(x -> sum(body_to_wind_axes(Line(x[:,1], x[:,2]), α, β).r1), x)

## ReverseDiff
ReverseDiff.gradient(x -> sum(body_to_wind_axes(x - r2, α, β)), r1) # (FAILS)
ReverseDiff.gradient(x -> let l = Line(x[:,1], x[:,2]); sum([ l.r1; l.r2 ]) end, x)
ReverseDiff.gradient(x -> let l = Line(x[:,1], x[:,2]); sum(body_to_wind_axes(sum([ l.r1 l.r2 ], dims = 2), α, β)) end, x)
ReverseDiff.gradient(x -> sum(body_to_wind_axes(Line(x[:,1], x[:,2]), α, β).r1 - r2), x)

## Zygote
Zygote.gradient(x -> sum(body_to_wind_axes(x - r2, α, β)), r1)
Zygote.gradient(x -> let l = Line(x[:,1], x[:,2]); sum([ l.r1; l.r2 ]) end, x)
Zygote.gradient(x -> let l = Line(x[:,1], x[:,2]); sum(body_to_wind_axes(sum([ l.r1 l.r2 ], dims = 2), α, β)) end, x)
Zygote.gradient(x -> sum(body_to_wind_axes(Line(x[:,1], x[:,2]), α, β).r1), x)