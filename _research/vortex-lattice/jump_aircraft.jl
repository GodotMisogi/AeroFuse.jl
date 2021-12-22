##
using JuMP
using Ipopt
using AeroMDAO
using ComponentArrays
using StaticArrays

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

function run_case(wing, htail_horses, vtail_horses, V, α, ρ, span_num, chord_num)
    uniform = Freestream(V, α, 0.0, zeros(3))

    wing_mesh  = WingMesh(wing, [span_num], chord_num, 
                          span_spacing = Cosine()
                         )
    
    aircraft = ComponentArray(
                                wing  = make_horseshoes(wing_mesh),
                                htail = htail_horses,
                                vtail = vtail_horses
                            );
    
    system  = solve_case(aircraft, uniform, refs)

    nearfield(system), farfield(system)
end

# Objective function
evaluate_CDi(wing, htail_horses, vtail_horses, V, α, ρ, span_num, chord_num) = run_case(wing, htail_horses, vtail_horses, V, α, ρ, span_num, chord_num)[1][1]

# Lift and area constraint function
function evaluate_cons(wing, htail_horses, vtail_horses, V, α, ρ, span_num, chord_num)
    area   = projected_area(wing)
    nf, ff = run_case(wing, htail_horses, vtail_horses, V, α, ρ, span_num, chord_num)
    lift   = dynamic_pressure(ρ, V) * area * nf[3]

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

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4(0,0,1,2), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4(0,0,0,9), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

htail_mesh = WingMesh(htail, [6], 4, 
                      span_spacing = Cosine()
                      )
vtail_mesh = WingMesh(vtail, [6], 4, 
                      span_spacing = Cosine()
                      )

htail_horses = make_horseshoes(htail_mesh)
vtail_horses = make_horseshoes(vtail_mesh)

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
ρ       = 1.225
ref     = SVector(x_w, 0., 0.)
refs    = References(S, b, c, ρ, ref)



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
evaluate_CDi(x...)  = evaluate_CDi(make_wing(x[1:end-1]), htail_horses, vtail_horses, V, x[end], ρ, span_num, chord_num)
evaluate_cons(x...) = evaluate_cons(make_wing(x[1:end-1]), htail_horses, vtail_horses, V, x[end], ρ, span_num, chord_num)
run_case(x...)      = run_case(make_wing(x[1:end-1]), htail_horses, vtail_horses, V, x[end], ρ, span_num, chord_num)

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