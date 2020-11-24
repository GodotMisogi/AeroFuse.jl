## 
using Revise
using JuMP
using Ipopt
using StaticArrays

##
includet("../src/FoilParametrization.jl")
using .FoilParametrization: naca4
using AeroMDAO

## Wing setup
#============================================#

# Parameters
α_u = [0.2, 0.3, 0.2, 0.15, 0.2]
dzs = (0., 0.)

foil = naca4((2,4,1,2))
num_secs = 3
foils = [ foil for i ∈ 1:num_secs ]

airfoils = Foil.(foils)
wing_twists = [2., 0., -2.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0., 11.3]
wing_sweeps = [1.14, 8.]

## Case setup
ρ = 1.225
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 0.0, 0.0)


function make_wing(chords...) 
    wing_right = HalfWing(airfoils, [ chords... ], wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
    Wing(wing_right, wing_right)
end

# Objective function
function optimize_cdi_chords(chords...)
    wing = make_wing(chords...)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    cdi = solve_case(wing, uniform, Ω, ref, span_num = 10, chord_num = 5, print = false)[2]
end

# Lift constraint function
function lift(chords...)
    wing = make_wing(chords...)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    cl = solve_case(wing, uniform, Ω, ref, span_num = 10, chord_num = 5, print = false)[1]
    lift = 0.5 * ρ * uniform.mag^2 * projected_area(wing) * cl
end

## Test runs
wing_chords = [0.18, 0.16, 0.08]
cdi = optimize_cdi_chords(wing_chords...)
lifter = lift(wing_chords...)
println("CDi: $cdi")
println("Lift: $lifter N")

## Optimization
wing_design = Model(Ipopt.Optimizer)
num_dv = 3

@variable(wing_design, 1e-4 <= chords[1:num_dv] <= 1.) 

register(wing_design, :optimize_chords, num_dv, optimize_cdi_chords, autodiff = true)

# register(wing_design, :lift_constraint, num_dv, lift, autodiff = true)

@NLobjective(wing_design, Min, optimize_chords(chords...))
# @NLconstraint(wing_design, lift_constraint(chords...) == 20)

## Run optimization
optimize!(wing_design)

## Print
println("Chords: $(value(chords)), Optimal CDi: $(objective_value(wing_design))")
chord_vals = value.(chords)