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
                         viscous = true)[2]
end

# Objective function
evaluate_CDi(wing, V, α, ρ, span_num, chord_num) = run_case(wing, V, α, ρ, span_num, chord_num)[1]

# Lift and area constraint function
function evaluate_cons(wing, V, α, ρ, span_num, chord_num)
    area = projected_area(wing)
    ff   = run_case(wing, V, α, ρ, span_num, chord_num)
    lift = dynamic_pressure(ρ, V) * area * ff[5]

    lift
end

## Test runs
#============================================#


# Parameters
V, α, ρ  = 29., 5., 1.225

# Design variables
n    = 5
global wing = Wing(foils     = fill(Foil(naca4(2,4,1,2)), n),
                   chords    = fill(0.314, n),
                   twists    = fill(0.0, n),
                   spans     = fill(1.3/(n-1), n-1),
                   dihedrals = fill(0., n-1),
                   LE_sweeps = fill(0., n-1))

# Meshing and assembly
wing_mac = mean_aerodynamic_center(wing);
b, S, c  = info(wing)[1:end-1]

span_num  = 10
chord_num = 1

x0 = [(chords ∘ right)(wing); α]

make_wing(x) = Wing(chords    = collect(x[1:5]),
                    foils     = foils(right(wing)),
                    spans     = spans(right(wing)),
                    twists    = rad2deg.(twists(right(wing))),
                    dihedrals = rad2deg.(dihedrals(right(wing))),
                    LE_sweeps = collect(x[6:end]))

# Closures
evaluate_CDi(x...)  = evaluate_CDi(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
evaluate_cons(x...) = evaluate_cons(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
run_case(x...)      = run_case(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)

##
wing_chords = [0.18, 0.16, 0.08, 0.04, 0.02]
wing_sweeps = [0, 5, 10, 20.]
x0 = ComponentVector(
                     chords = wing_chords,
                     sweeps = wing_sweeps,
                     alpha  = α
                     )
wing = make_wing(x0)

##
cdi    = evaluate_CDi(x0...)
# lifter = evaluate_cons(x0...)
# coeffs = run_case(x0...)
# println("CDi: $cdi")
# println("CL: $(force_coefficient(lifter, dynamic_pressure(1.225, 10.), projected_area(wing)))")

## Chord optimization
#============================================#

chord_design = Model(with_optimizer(Ipopt.Optimizer))
num_dv = 10

@variable(chord_design, chordsweepsalpha[i=1:num_dv], start = x0[i])

# set_value(chordsweepsalpha, x0)

##
register(chord_design, :evaluate_CDi, num_dv, evaluate_CDi, autodiff=true)
# register(chord_design, :evaluate_cons, num_dv, evaluate_cons, autodiff = true)

##
@NLobjective(chord_design, Min, evaluate_CDi(chordsweepsalpha...))
# @NLconstraint(chord_design, evaluate_cons(chordsweepsalpha...) == 10.)

## Run optimization
optimize!(chord_design)

## Print optimal case
println("Chords: $(value.(chordsweepsalpha[1:5])), Optimal CDi: $(objective_value(chord_design))")
chord_vals = value.(chords)
uniform = Freestream(10., 0., 0.)
wing = make_wing(chord_vals...)
nf_coeffs, ff_coeffs, horseshoe_panels, normals, horseshoes, Γs = run_case(wing, 10., 0., 0.)
cdi = ff_coeffs[2]
cl = ff_coeffs[1]
begin
    println("\nNearfield:")
    print_dynamics(nf_coeffs...)
    println("\nFarfield:")
    print_dynamics(ff_coeffs...)
end

##
layout = Layout(title = "Vortex Lattice",
                scene = attr(aspectratio=attr(x=1,y=1,z=1)));

## Streamlines
num_points = 30
max_z = 0.02
y = span(wing) / 2 - 0.5
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

span_points = 10
init = trailing_chopper(wing.right, span_points)
dx, dy, dz = 0, 0, 1e-3
seed = [ init .+ Ref([dx, dy, dz]);
         init .+ Ref([dx, dy, -dz]) ]

trace_streams = trace_streamlines(uniform, seed, horseshoes[:], Γs[:], 2, 100);

##
trace_horsies = trace_panels(horseshoe_panels[:])
trace_horses = trace_panels(horseshoe_panels[:], Γs[:])
trace_wing = trace_surface(wing)

PlotlyJS.plot(
            [
                (trace for trace in trace_horsies)...,
                (trace for trace in trace_horses)...,
                (trace for trace in trace_streams)...,
                # [ trace for trace in trace_wing ]...,
            ],
            layout)

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