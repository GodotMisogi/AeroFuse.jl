## 
using Revise
using StaticArrays
using AeroMDAO
using Optim
using PlotlyJS

# Wing setup
#============================================#

## Helper functions
function make_wing(digits, chords, twists, spans, dihedrals, sweeps)
    wing_right = HalfWing(Foil.(naca4(tuple(digits...)) for i ∈ 1:5), chords, twists, spans, dihedrals, sweeps)
    Wing(wing_right, wing_right)
end


function run_case(wing, V, α)
    uniform = Freestream(V, α, 0.0)
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 6, print = false)

    coeffs
end

# Objective function
optimize_cdi(digits, chords, twists, spans, dihedrals, sweeps, V, α) = run_case(make_wing(digits, chords, twists, spans, dihedrals, sweeps), V, α)[2]

# Lift constraint function
cons_lift(digits, chords, twists, spans, dihedrals, sweeps, V, α, ρ = 1.225) = 
    let wing = make_wing(digits, chords, twists, spans, dihedrals, sweeps); dynamic_pressure(ρ, V) * projected_area(wing) * run_case(wing, V, α)[1] end

## Test runs
#============================================#

# Parameters
digits = [2,4,1,2]
chords = [0.18, 0.16, 0.08, 0.04, 0.02]
twists = [0., 0., 0., 0., 0.]
spans = [0.2, 0.2, 0.2, 0.2]
dihedrals = [0., 0., 0., 0.]
sweeps = [0., 0., 0., 0.]
V = 10.
α = 0.

wing = make_wing(digits, chords, twists, spans, dihedrals, sweeps)
cdi = optimize_cdi(digits, chords, twists, spans, dihedrals, sweeps, V, α)
lifter = cons_lift(digits, chords, twists, spans, dihedrals, sweeps, V, α)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, V), projected_area(wing)))")

## Optimisation
#============================================#

optimize_chords = TwiceDifferentiable(x -> optimize_cdi(digits, x, twists, spans, dihedrals, sweeps, V, α), chords, autodiff = :forward)

l_bound = fill(1e-12, length(chords))
u_bound = fill(Inf, length(chords))
bound_chords = TwiceDifferentiableConstraints(l_bound, u_bound)

lc = [5.0]; uc = [5.0]
lift_constraint = TwiceDifferentiableConstraints(x -> cons_lift(digits, x, twists, spans, dihedrals, sweeps, V, α), l_bound, u_bound, lc, uc, :forward)

resi_chords = optimize(optimize_chords,
                        lift_constraint,
                        chords,							# Initial value
                        IPNewton(),
                        autodiff = :forward,
                        Optim.Options(
                                        # extended_trace = true,
                                        show_trace = true
                                     )
                        )