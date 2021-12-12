##
using JuMP
using Ipopt
using AeroMDAO
using ComponentArrays
# using PlotlyJS

# Wing setup
#============================================#

## Helper functions
function make_wing(x)
    wing_right = HalfWing(Foil.(naca4((0,0,1,2)) for i = 1:5),
                          [ x[1:5]... ],
                          [0., 0., 0., 0., 0.],
                          [0.2, 0.2, 0.2, 0.2],
                          [0., 0., 0., 0.],
                          [ x[6:end]... ])
    Wing(wing_right, wing_right)
end

## FUNCTIONAL TEST
#==========================================================================================#

# Helper functions
#============================================#

function run_case(wing, V, α, ρ, span_num, chord_num)
    uniform = Freestream(V, α, 0.0, zeros(3))
    coeffs  = solve_case(wing, uniform,
                         rho_ref = ρ,
                         span_num = span_num,
                         chord_num = chord_num,
                         viscous = true)[1:2]
end

# Objective function
evaluate_CDi(wing, V, α, ρ, span_num, chord_num) = run_case(wing, V, α, ρ, span_num, chord_num)[1][1]

# Lift and area constraint function
function evaluate_cons(wing, V, α, ρ, span_num, chord_num)
    area   = projected_area(wing)
    nf, ff = run_case(wing, V, α, ρ, span_num, chord_num)
    lift   = dynamic_pressure(ρ, V) * area * nf[5]

    lift
end

## Test runs
#============================================#


# Parameters
V, α, ρ  = 15., 5., 1.225

# Design variables
n    = 12
wing = Wing(foils     = fill(Foil(naca4(2,4,1,2)), n),
            chords    = fill(0.314, n),
            twists    = fill(0.0, n),
            spans     = fill(1.3/(n-1), n-1),
            dihedrals = fill(0., n-1),
            LE_sweeps = fill(0., n-1))

# Meshing and assembly
wing_mac = mean_aerodynamic_center(wing);
b, S, c  = info(wing)[1:end-1]

span_num  = 12
chord_num = 6

x0 = [chords(wing.right); sweeps(wing.right); α]

make_wing(x) = Wing(chords    = collect(x[1:n]),
                    foils     = foils(wing.right),
                    spans     = spans(wing.right),
                    twists    = rad2deg.(twists(wing.right)),
                    dihedrals = rad2deg.(dihedrals(wing.right)),
                    LE_sweeps = collect(x[n+1:end]))

# Closures
evaluate_CDi(x...)  = evaluate_CDi(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
evaluate_cons(x...) = evaluate_cons(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
run_case(x...)      = run_case(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)

##
wing_chords = LinRange(0.2, 0.1, n)
wing_sweeps = fill(10, n-1)
x0 = ComponentVector(
                     chords = wing_chords,
                     sweeps = wing_sweeps,
                     alpha  = α
                     )
wing = make_wing(x0)

## Test runs
@show cdi    = evaluate_CDi(x0...)
@show lifter = evaluate_cons(x0...)
@show coeffs = run_case(x0...)

## Chord optimization
#============================================#

chord_design = Model(with_optimizer(Ipopt.Optimizer))
set_optimizer_attribute(chord_design, "max_iter", 100)

## Variables and bounds
num_dv = 2n
@variable(chord_design, 0. <= x[i=1:num_dv] <= 20, start = x0[i])

## Registration
register(chord_design, :evaluate_CDi, num_dv, evaluate_CDi, autodiff=true)
register(chord_design, :evaluate_cons, num_dv, evaluate_cons, autodiff = true)

## Objective
@NLobjective(chord_design, Min, evaluate_CDi(x...))

## Constraints
@NLconstraint(chord_design, evaluate_cons(x...) == 35.)

## Run optimization
optimize!(chord_design)

## Print optimal case
@show value.(x)
@show objective_value(chord_design)

##
x_opt = value.(x)
wing = make_wing(x_opt[1:end-1]);

##
nf, ff = run_case(wing, V, x_opt[end], ρ, span_num, chord_num)
print_coefficients(nf, ff)

## Plotting
bing = plot_wing(wing)

##
using Plots

plot(bing)

## NACA-4 optimization
#============================================#

naca4_design = Model(with_optimizer(Ipopt.Optimizer))

@variable(naca4_design, 1e-4 <= digits[1:4] <= 1.)

register(naca4_design, :optimize_cdi_naca4, 4, optimize_cdi_naca4, autodiff = true)

register(naca4_design, :naca4_lift, 4, naca4_lift, autodiff = true)

@NLobjective(naca4_design, Min, optimize_cdi_naca4(digits...))
@NLconstraint(naca4_design, naca4_lift(digits...) == 5)

## Run optimization
optimize!(naca4_design)

## Print
println("Digits: $(value.(chords)), Optimal CDi: $(objective_value(naca4_design))")
digit_vals = value.(digits)
cdi = optimize_cdi_chords(digit_vals...)
lifter = lift(digit_vals...)
println("CDi: $cdi")
println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")