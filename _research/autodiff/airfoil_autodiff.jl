using AeroMDAO
using BenchmarkTools

##
using ForwardDiff, ReverseDiff, Zygote

using StaticArrays
# @Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
# @Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)

## Geometry tests
n_dim = 200

xs = naca4((2,4,1,2))
xi = cosine_foil(xs, n_dim)

foil = Foil(xi[:,1], xi[:,2])
len  = arc_length(foil)

##
@time ForwardDiff.gradient((arc_length ∘ Foil), xi)

##
@time ReverseDiff.gradient((arc_length ∘ Foil), xi)

##
Zygote.gradient(arc_length ∘ Foil, xi)
Zygote.gradient(arc_length ∘ Foil, xi[:,1], xi[:,2])

## Aerodynamic analysis
function optimize_CST(coords)
    # airfoil = Foil(coords)
    # airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), (0., 0.), 80)
    vecs   = SVector.(coords[:,1], coords[:,2])
    panels = Panel2D.(vecs[2:end], vecs[1:end-1])[end:-1:1]
    AeroMDAO.solve_problem(panels, sincos(deg2rad(5.))[end:-1:1], deg2rad(5.), false, 1e5)[end]
end

alpha_u = fill(0.1, 60)
alpha_l = fill(-0.1, 60)

foil = Foil(xi)

optimize_CST(xi)

## PASSES
@time ForwardDiff.gradient(optimize_CST, xi)

## PASSES
@time ReverseDiff.gradient(optimize_CST, xi)

## FAILS (due to binomial and SVector)
# @time Zygote.jacobian(optimize_CST, xs)