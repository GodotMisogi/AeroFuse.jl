using ForwardDiff, ReverseDiff, Zygote
using LinearAlgebra

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