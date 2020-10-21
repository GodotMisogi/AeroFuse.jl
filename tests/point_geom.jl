##
using Revise
includet("../src/MathTools.jl")
includet("../src/Geometry.jl")
includet("../src/LiftingLine.jl")
includet("../src/AeroMDAO.jl")
includet("../src/FoilParametrization.jl")

##
import Base: *, +
using LinearAlgebra
using StaticArrays
using .AeroMDAO: Foil
using .LiftingLine: Panel3D
using .FoilParametrization: read_foil

function svectors(x::Matrix{T}, ::Val{N}) where {T,N}
    size(x,1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    reinterpret(SVector{N,T}, vec(x))
end

## Wing section setup
foilpath = "airfoil_database/ys930.dat"

num_secs = 6

coords = read_foil(foilpath)
# foils = [ coords for n in 1:num_secs ]
airfoil = Foil(coords)
2 * airfoil.coords

##
xs = [1 2 3; 4 5 6; 7 8 9]
ys = [1 0 1; 1 1 0; 0 1 1]' * xs
zs = SArray{3}(xs)


##
pans = [ ys[1:end-1,:] ys[2:end,:] xs[2:end,:] xs[1:end-1,:] ]

##
[ Panel3D(ps...) for ps in eachrow(pans) ]

##
