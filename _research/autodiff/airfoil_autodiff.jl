using AeroMDAO
using BenchmarkTools

##
using ForwardDiff, ReverseDiff, Zygote

# using StaticArrays
# @Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
# @Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)

## Geometry tests
n_dim = 200

xs = naca4((2,4,1,2))
xi = cosine_foil(xs, n_dim)

foil = Foil(xi[:,1], xi[:,2])
len  = arc_length(foil)

##
@benchmark ForwardDiff.gradient($(arc_length ∘ Foil), $xi)

##
@benchmark ReverseDiff.gradient($(arc_length ∘ Foil), $xi)

##
Zygote.gradient(arc_length ∘ Foil, xi)
Zygote.gradient(arc_length ∘ Foil, xi[:,1], xi[:,2])

## Aerodynamic analysis
function optimize_CST(coords)
    airfoil = Foil(coords)
    # airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), (0., 0.), 80)
    uniform = Uniform2D(1.0, 5.0)
    solve_case(airfoil, uniform, num_panels = size(coords)[1])[1]
end

alpha_u = fill(0.1, 60)
alpha_l = fill(-0.1, 60)

foil = Foil(xi)

optimize_CST(xi)

## PASSES
@benchmark ReverseDiff.gradient($optimize_CST, $xi)

## PASSES
@benchmark ForwardDiff.gradient($optimize_CST, $xi)

## FAILS (due to binomial and SVector)
# @time Zygote.jacobian(optimize_CST, xs)