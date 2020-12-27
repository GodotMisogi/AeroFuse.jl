## 
using Revise
using JuMP
using Ipopt
using StaticArrays
using AeroMDAO

# Wing setup
#============================================#

## Helper functions
function make_wing(chords...)
    wing_right = HalfWing(Foil.(naca4((2,4,1,2)) for i ∈ 1:5), 
                         [ chords... ], 
                         [0., 0., 0., 0., 0.],
                         [0.2, 0.2, 0.2, 0.2],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 0.])
    Wing(wing_right, wing_right)
end

function naca4_wing(m, p, t, c)
    wing_right = HalfWing(Foil.(naca4((m, p, t, c)) for i ∈ 1:5), 
                         [0.18, 0.16, 0.08, 0.04, 0.02],
                         [0., 0., 0., 0., 0.],
                         [0.2, 0.2, 0.2, 0.2],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 0.])
    Wing(wing_right, wing_right)
end

# Objective function
function optimize_cdi_chords(chords...)
    wing = make_wing(chords...)
    ρ = 1.225
    Ω = SVector(0.0, 0.0, 0.0)
    uniform = Freestream(10.0, 0.0, 0.0, Ω)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, cps, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10)

    cdi = ff_coeffs[2]
end

function optimize_cdi_naca4(m, p, t, c)
    wing = naca4_wing(m, p, t, c)
    ρ = 1.225
    Ω = SVector(0.0, 0.0, 0.0)
    uniform = Freestream(10.0, 0.0, 0.0, Ω)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, cps, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10)

    cdi = ff_coeffs[2]
end

# Lift constraint function
function chords_lift(chords...)
    wing = make_wing(chords...)
    ρ = 1.225
    Ω = SVector(0.0, 0.0, 0.0)
    uniform = Freestream(10.0, 0.0, 0.0, Ω)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, cps, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10)

    cl = ff_coeffs[1]
    lift = 1/2 * ρ * uniform.mag^2 * projected_area(wing) * cl
end

function naca4_lift(m, p, t, c)
    wing = naca4_wing(m, p, t, c)
    ρ = 1.225
    Ω = SVector(0.0, 0.0, 0.0)
    uniform = Freestream(10.0, 0.0, 0.0, Ω)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, cps, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10)

    cl = ff_coeffs[1]
    lift = 1/2 * ρ * uniform.mag^2 * projected_area(wing) * cl
end

## Test runs
#============================================#

wing_chords = [0.18, 0.16, 0.08, 0.04, 0.02]
wing = make_wing(wing_chords...)
cdi = optimize_cdi_chords(wing_chords...)
lifter = chords_lift(wing_chords...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

##
naca_test = (2,4,1,2)
wing = naca4_wing(naca_test...)
cdi = optimize_cdi_naca4(naca_test...)
lifter = naca4_lift(naca_test...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

## Chord optimization
#============================================#

chord_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 5

@variable(chord_design, 1e-4 <= chords[1:num_dv] <= 1.) 

register(chord_design, :optimize_cdi_chords, num_dv, optimize_cdi_chords, autodiff = true)

register(chord_design, :chords_lift, num_dv, chords_lift, autodiff = true)

@NLobjective(chord_design, Min, optimize_cdi_chords(chords...))
@NLconstraint(chord_design, chords_lift(chords...) == 5)

## Run optimization
optimize!(chord_design)

## Print
println("Chords: $(value.(chords)), Optimal CDi: $(objective_value(chord_design))")
chord_vals = value.(chords)
cdi = optimize_cdi_chords(chord_vals...)
lifter = chords_lift(chord_vals...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

## NACA-4 optimization
#============================================#

naca4_design = Model(with_optimizer(Ipopt.Optimizer))

@variable(naca4_design, 1e-4 <= digits[1:4] <= 1.) 

register(naca4_design, :optimize_cdi_naca4, 4, optimize_cdi_naca4, autodiff = true)

register(naca4_design, :naca4_lift, 4, naca4_lift, autodiff = true)

@NLobjective(naca4_design, Min, optimize_cdi_naca4(digits...))
@NLconstraint(naca4_design, naca4_lift(digits...) == 0.6)

## Run optimization
optimize!(naca4_design)

## Print
println("Digits: $(value.(chords)), Optimal CDi: $(objective_value(naca4_design))")
digit_vals = value.(digits)
cdi = optimize_cdi_chords(digit_vals...)
lifter = lift(digit_vals...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")