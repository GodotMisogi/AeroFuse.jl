using AeroMDAO
using BenchmarkTools

##
using ForwardDiff, Zygote
using ReverseDiff: GradientTape, GradientConfig, compile, DiffResults, gradient, gradient!
using StaticArrays
# @Zygote.adjoint (T :: Type{<:SVector})(xs :: Number ...) = T(xs...), dv -> (nothing, dv...)
# @Zygote.adjoint (T :: Type{<:SVector})(x :: AbstractVector) = T(x), dv -> (nothing, dv)

## Geometry tests
n_dim = 200

xs = naca4((2,4,1,2))
xi = cosine_spacing(xs, n_dim)

foil = Foil(xi[:,1], xi[:,2])
len  = arc_length(foil)

##
@time ForwardDiff.gradient((arc_length ∘ Foil), xi)

##
@time gradient((arc_length ∘ Foil), xi)

##
Zygote.gradient(arc_length ∘ Foil, xi)
Zygote.gradient(arc_length ∘ Foil, xi[:,1], xi[:,2])

## Aerodynamic analysis
function optimize_CST(coords)
    airfoil = Foil(coords)
    fs = Uniform2D(1., 5.)
    solve_case(airfoil, fs, num_panels = length(airfoil.x) - 1)[1]
    # airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), (0., 0.), 80)
    # vecs   = SVector.(coords[:,1], coords[:,2])
    # panels = Panel2D(coords)
    
    # AeroMDAO.solve_problem(panels, sincos(deg2rad(5.))[end:-1:1], deg2rad(5.), false, 1e5)[end]
end

alpha_u = fill(0.1, 60)
alpha_l = fill(-0.1, 60)

foil = Foil(xi)

##
@time optimize_CST(xi)

## PASSES
@time ForwardDiff.gradient(optimize_CST, xi)

## PASSES
@time gradient(optimize_CST, xi)

## FAILS (due to binomial and SVector)
# @time Zygote.jacobian(optimize_CST, xs)

## Taping
const f_tape = GradientTape(optimize_CST, xi)

# compile `f_tape` into a more optimized representation
const compiled_f_tape = compile(f_tape)

##
results = similar(xi)
all_results = DiffResults.GradientResult(results)
cfg = GradientConfig(xi)

##
@time gradient!(results, optimize_CST, xi, cfg)

##
@time gradient!(all_results, optimize_CST, xi, cfg)

##
@time gradient!(results, compiled_f_tape, xi)