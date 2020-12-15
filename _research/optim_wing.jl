## 
using Revise
using JuMP
using Ipopt
using StaticArrays
using AeroMDAO

# Wing setup
#============================================#

## Helper functions
function make_wing(chords)
    wing_right = HalfWing(Foil.(naca4((2,4,1,2)) for i ∈ 1:5), 
                         chords, 
                         [0., 0., 0., 0., 0.],
                         [0.2, 0.2, 0.2, 0.2],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 0.])
    Wing(wing_right, wing_right)
end

function naca4_wing(digits)
    wing_right = HalfWing(Foil.(naca4(tuple(digits)) for i ∈ 1:5), 
                         [0.18, 0.16, 0.08, 0.04, 0.02],
                         [0., 0., 0., 0., 0.],
                         [0.2, 0.2, 0.2, 0.2],
                         [0., 0., 0., 0.],
                         [0., 0., 0., 0.])
    Wing(wing_right, wing_right)
end

function run_chords(chords, V, α)
    wing = make_wing(chords)
    uniform = Freestream(V, α, 0.0)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 6, print = false)

    coeffs
end

# Objective function
optimize_cdi_chords(chords, V, α) = run_chords(chords, V, α)[2]

# Lift constraint function
chords_lift(chords, V, α, S, ρ = 1.225) = dynamic_pressure(ρ, V) * S * run_chords(chords, V, α)[1]

## Test runs
#============================================#

wing_chords = [0.18, 0.16, 0.08, 0.04, 0.02]
V = 10.
α = 0.
wing = make_wing(wing_chords)
cdi = optimize_cdi_chords(wing_chords, V, α)
lifter = chords_lift(wing_chords, V, α, projected_area(wing))
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, V), projected_area(wing)))")

##
naca_test = (2,4,1,2)
wing = naca4_wing(naca_test)
cdi = optimize_cdi_naca4(naca_test)
lifter = naca4_lift(naca_test)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, V), projected_area(wing)))")