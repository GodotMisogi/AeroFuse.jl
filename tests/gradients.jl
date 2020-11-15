## 
using Revise
includet("../src/FoilParametrization.jl")

##
using .FoilParametrization: naca4, kulfan_CST
using AeroMDAO
using BenchmarkTools
using ProfileView
using StaticArrays

# Gradient package
using Zygote


## 2D doublet-source panel method
alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.2, -0.1, -0.1, -0.001]
dzs = (1e-4, 1e-4)

uniform = Uniform2D(1.0, 5.0)
airfoil = kulfan_CST(alpha_u, alpha_l, (0., 0.), 0.0, 80)
panels = make_2Dpanels(airfoil)


function optimize_CST(α)
    cl = solve_case(panels, Uniform2D(1.0, α))
end

optimize_CST(5.0)

##
@Zygote.adjoint (T::Type{<:SVector})(xs::Number...) = T(xs...), dv -> (nothing, dv...)
@Zygote.adjoint (T::Type{<:SVector})(x::AbstractVector) = T(x), dv -> (nothing, dv)

gradient(optimize_CST, 3.0)


## 3D vortex-lattice method
num_secs = 3

foil = naca4((4,4,1,2))
foils = [ foil for i ∈ 1:num_secs ]
airfoils = Foil.(foils)
wing_twists = [2., 0., -2.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0, 11.3]
wing_sweeps = [1.14, 8]

ρ = 1.225
uniform = Uniform(10.0, 5.0, 0.0)

## Define objective
function test_model(wing_chords)
    wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
    wing = Wing(wing_right, wing_right)

    ## Assembly
    ref = (0.25 * mean_aerodynamic_chord(wing), 0, 0)
    coeffs = solve_case(wing, uniform, ref, span_num = 10, chord_num = 5, print = false)

    return coeffs
end

lift_coefficient(chords) = test_model(chords)[1]
drag_coefficient(chords) = test_model(chords)[2]
lift_to_drag_ratio(chords) = lift_coefficient(chords)/drag_coefficient(chords)

##
lift_coefficient([0.16, 0.12, 0.08])

##
lift_coefficient([0.2, 0.1, 0.05])


##
using Zygote

lift_coefficient'([0.16, 0.12, 0.08])
# grad_lift = inputs -> gradient(lift_coefficient, inputs)

## Gradient testing
# using ReverseDiff: GradientTape, GradientConfig, gradient, gradient!, compile, DiffResults

# ## pre-record a GradientTape for `f` using inputs of shape 3 with Float64 elements
# const cd_tape = GradientTape(drag_coefficient, rand(3))

# ## compile `f_tape` into a more optimized representation
# const compiled_cd_tape = compile(cd_tape)

# ##
# x = rand(3)
# inputs = (x)
# results = similar(x)
# all_results = DiffResults.GradientResult(results)
# cfg = GradientConfig(inputs)

## Optimization libraries
# using JuMP
# using Ipopt

# wing_design = Model(Ipopt.Optimizer)

# ##
# @variable(wing_design, 1e-4 <= root_chord <= 0.5, start = 0.16)
# @variable(wing_design, 1e-4 <= rootip_chord <= 0.5, start = 0.12)
# @variable(wing_design, 1e-4 <= tip_chord <= 0.5, start = 0.08)

# register(wing_design, :test_model, 1, drag_coefficient([root_chord, rootip_chord, tip_chord]), autodiff = true)
# register(wing_design, :lift_coefficient, 1, lift_coefficient([root_chord, rootip_chord, tip_chord]), autodiff = true) 

# @NLconstraint(wing_design, con, lift_coefficient([root_chord, rootip_chord, tip_chord]) == 0.8)

# @NLobjective(wing_design, Min, drag_coefficient([root_chord, rootip_chord, tip_chord]))

# ##
# optimize!(wing_design)

# ##
# objective_value(wing_design)
# vals = value.([root_chord, rootip_chord, tip_chord])
