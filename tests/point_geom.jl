##
include("../src/MathTools.jl")
include("../src/Geometry.jl")
include("../src/LiftingLine.jl")

##
import Base: *, +
using .LiftingLine: Panel3D
using LinearAlgebra
using StaticArrays

##
pts = [1.0 2 3; 4 5 6; 7 8 9]
xs = SVector.(pts[:,1], pts[:,2], pts[:,3])

##
ys = [1 0 1; 1 1 0; 0 1 1]' * xs

##
pans = [ ys[1:end-1] ys[2:end] xs[2:end] xs[1:end-1] ]

##
[ Panel3D(ps...) for ps in eachrow(pans) ]