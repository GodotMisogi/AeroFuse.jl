## GalacticOptim.jl tests
using AeroMDAO
using GalacticOptim, Optim

# Helper functions
#============================================#

function run_case(wing, α, V, ρ, span_num, chord_num)
    fs = Freestream(alpha = α, beta = 0.0, omega = zeros(3))
    refs = References(speed = V, density = ρ, area = projected_area(wing), span = span(wing), chord = mean_aerodynamic_chord(wing), location = mean_aerodynamic_center(wing))
    wing_mesh = WingMesh(wing, fill(span_num ÷ length(wing.right.spans), length(wing.right.spans)), chord_num)
    aircraft = ComponentVector(wing = make_horseshoes(wing_mesh));
    system = solve_case(aircraft, fs, refs;)
end

# Objective function
evaluate_CDi(wing, α, V, ρ, span_num, chord_num) = farfield(run_case(wing, α, V, ρ, span_num, chord_num))[1]

# Lift and area constraint function
function evaluate_cons(wing, α, V, ρ, span_num, chord_num)
    sys  = run_case(wing, α, V, ρ, span_num, chord_num)
    ff   = farfield(sys)
    lift = dynamic_pressure(sys.reference.density, sys.reference.speed) * sys.reference.area * ff[3]
    # CDi  = ff[1]

    [lift, sys.reference.area]
end

## Optimization
#============================================#

# Parameters
V, α, ρ  = 29., 5., 1.225

# Design variables
n    = 4
wing = Wing(foils     = fill(naca4((2,4,1,2)), n),
            chords    = fill(0.314, n),
            twists    = fill(0.0, n),
            spans     = fill(1.3/(n-1), n-1),
            dihedrals = fill(0., n-1),
            sweeps      = fill(0., n-1))

make_wing(x) = Wing(chords    = x,
                    foils     = foils(wing.right),
                    spans     = spans(wing.right),
                    twists    = rad2deg.(twists(wing.right)),
                    dihedrals = rad2deg.(dihedrals(wing.right)),
                    sweeps      = rad2deg.(sweeps(wing.right)))

# Meshing and assembly
AR, b, S, c, mac  = properties(wing)

span_num  = 10
chord_num = 1

## Bounds and constraints
weight      = 12 * 9.81
load_factor = 1.5
lift_req    = weight * load_factor

l_bound = fill(1e-12, length(x0))
u_bound = fill(Inf,   length(x0))

lc = [ lift_req, projected_area(wing), 0. ]
uc = [ lift_req, projected_area(wing), Inf ]

# Initial guess and parameters
x0 = [ wing.right.chords; α]
p  = (V, ρ, span_num, chord_num)


# Closures and test
evaluate_CDi(x, p) = evaluate_CDi(make_wing(x[1:end-1]), x[end], p...)
cons(x, p) = evaluate_cons(make_wing(x[1:end-1]), x[end], p...)

test_CDi = @time evaluate_CDi(x0, p)
test_con = @time cons(x0, p)

## GalacticOptim variables
prob = OptimizationFunction(evaluate_CDi, GalacticOptim.AutoForwardDiff())
prob = OptimizationProblem(prob, x0, p) # lcons = l_bound, ucons = u_bound, in_place = false)

## Run optimisation
res_func = solve(prob, Newton())

## Sanity check
x_opt             = res_func.minimizer
CDi_opt, cons_opt = evaluate_CDi(x_opt), cons!(zeros(3), x_opt)

## Plotting
#============================================#

using Plots
pyplot(dpi = 300)

ini_plan = plot_wing(wing)

opt_wing = make_wing(x_opt[1:end-1])
opt_plan = plot_wing(opt_wing)

b = span(opt_wing)
plot(xlim = (-b/2, b/2), zlim = (-b/2, b/2))
plot!(ini_plan, label = "Initial")
plot!(opt_plan, color = :blue, label = "Optimized")


## Plotting
#============================================#

opt_state_wing = make_wing(x_opt[1:end-1])
opt_state_plan = plot_wing(opt_state_wing)
htail_plan     = plot_wing(htail,
                           position = [1., 0, 0.15],
                           angle    = deg2rad(0.),
                           axis     = [0., 1., 0.])
vtail_plan     = plot_wing(vtail_1,
                           position = [1., 0., 0],
                           angle    = π/2,
                           axis     = [1., 0, 0])

b = span(opt_wing)
plot(xlim = (-b/2, b/2), zlim = (-b/2, b/2))
plot!(ini_plan, label = "Initial")
plot!(opt_state_plan, color = :blue, label = "Optimized")
plot!(htail_plan, label = "HTail")
plot!(vtail_plan, label = "VTail")