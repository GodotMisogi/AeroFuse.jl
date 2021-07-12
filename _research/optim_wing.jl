## Optimisation tests
using Revise
using StaticArrays
using AeroMDAO
using Optim

## FUNCTIONAL TEST
#==========================================================================================#

# Helper functions
#============================================#

function run_case(wing, V, α, ρ, span_num, chord_num)
    uniform = Freestream(V, α, 0.0, zeros(3))
    coeffs  = solve_case(wing, uniform, 
                         rho_ref = ρ,
                         span_num = span_num,
                         chord_num = chord_num)[2]
end

# Objective function
evaluate_CDi(wing, V, α, ρ, span_num, chord_num) = run_case(wing, V, α, ρ, span_num, chord_num)[1]

# Lift and area constraint function
function evaluate_cons(wing, V, α, ρ, span_num, chord_num)
    area = projected_area(wing)
    lift = dynamic_pressure(ρ, V) * S * run_case(wing, V, α, ρ, span_num, chord_num)[3]

    lift, area
end

## Test runs
#============================================#

# Parameters
n    = 4
wing = Wing(foils     = fill(Foil(naca4((2,4,1,2))), n),
            chords    = fill(0.3, n),
            twists    = fill(0.0, n),
            spans     = fill(1.3/(n-1), n-1),
            dihedrals = fill(0.0, n-1),
            sweep_LEs = fill(10., n-1))

wing = WingSection(root_foil  = naca4((2,4,1,2)),
                   tip_foil   = naca4((2,4,1,2)),
                   span       = 2.6,
                   root_chord = 0.314,
                   taper      = 0.8,
                   root_twist = 0.0,
                   tip_twist  = 0.0,
                   dihedral   = 5.0,
                   sweep_LE   = 10.)

# Meshing and assembly
wing_mac = mean_aerodynamic_center(wing);
b, S, c  = info(wing)[1:end-1]
V, α, ρ  = 30., 4., 1.225

span_num  = 10
chord_num = 6

# Test case
ff_coeffs = @time evaluate_CDi(wing, V, α, ρ, span_num, chord_num)
cons_test = @time evaluate_cons(wing, V, α, ρ, span_num, chord_num)

## Optimisation
#============================================#

# Design variables - chords and angle of attack
x0           = [(chords ∘ right)(wing); α]

x0           = [mean_aerodynamic_chord(wing);
                taper_ratio(right(wing));
                dihedrals(right(wing));
                sweeps(right(wing));
                α]

weight       = 15 * 9.81
load_factor  = 1.5
lift_req     = weight * load_factor

# Bounds
l_bound      = fill(1e-12, length(x0))
u_bound      = fill(Inf,   length(x0))

# Constraint bounds
lc = [ lift_req, projected_area(wing) ]
uc = [ lift_req, projected_area(wing) ]

# Closures
make_wing(x) = Wing(chords    = x, 
                    foils     = foils(right(wing)),
                    spans     = spans(right(wing)),
                    twists    = rad2deg.(twists(right(wing))), 
                    dihedrals = rad2deg.(dihedrals(right(wing))),
                    sweep_LEs = rad2deg.(sweeps(right(wing))))

make_wing(x) = WingSection(root_foil  = naca4((2,4,1,2)),
                           tip_foil   = naca4((2,4,1,2)),
                           span       = span(wing),
                           root_chord = x[1],
                           taper      = x[2],
                           dihedral   = rad2deg(x[3]),
                           sweep_LE   = rad2deg(x[4]))

evaluate_CDi(x) = evaluate_CDi(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
cons!(cs, x) = cs .= evaluate_cons(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)

# Optim variables
optimize_chords = TwiceDifferentiable(evaluate_CDi, x0)
cons_lift_area  = TwiceDifferentiableConstraints(cons!, l_bound, u_bound, lc, uc);

## Run optimisation
res_chord = optimize(optimize_chords,          # Objective functions
                     cons_lift_area,           # Constraint
                     x0,                       # Initial value
                     IPNewton(),               # Optimization algorithm
                    #  autodiff = :forward,    # Automatic differentiation
                     Optim.Options(
                                #    extended_trace = true,
                                   show_trace = true
                                  )
                    )

## Plotting
#============================================#

using Plots
pyplot(dpi = 300)

ini_plan = plot_wing(wing)

opt_wing = make_wing(res_chord.minimizer[1:end-1])
opt_plan = plot_wing(opt_wing)

b = span(opt_wing)
plot(xlim = (-b/2, b/2), zlim = (-b/2, b/2))
plot!(ini_plan, label = "Initial")
plot!(opt_plan, color = :blue, label = "Optimized")

## STATEFUL TEST
#==========================================================================================#

# Helper functions
#============================================#

function make_surface(wing, span_num, chord_num)
    wing_panels, wing_normals = panel_wing(wing, span_num, chord_num;)
    VLMSurface(wing_panels, wing_normals, "Wing")
end

function run_case!(surf, α, system, state)
    # Update angle of attack and surface
    state.alpha        = α
    system.surfaces[1] = surf

    # Run case
    evaluate_case!(system, state)
end

# Objective function
evaluate_induced_drag(system) = sum(farfield_forces, surfaces(system))[1]

# Lift and area constraint function
evaluate_constraints(wing, system) = [ projected_area(wing); sum(farfield_forces, surfaces(system))[3] ]

## Optimisation
#============================================#

# Set up state
state = VLMState(V, deg2rad(α), 0., [0., 0., 0.], 
                 rho_ref   = ρ,
                 r_ref     = [ wing_mac[1], 0, 0 ],
                 area_ref  = projected_area(wing), 
                 chord_ref = mean_aerodynamic_chord(wing), 
                 span_ref  = span(wing));

# Solve initial case
aircraft = Dict("Wing" => panel_wing(wing, span_num, chord_num));
system   = solve_case(aircraft, state);

# Test case - Fixed speed
state.speed   = V
state.alpha   = deg2rad(α)
state.rho_ref = ρ

# Bounds
l_bound = fill(1e-12, length(x0))
u_bound = fill(Inf, length(x0))

# Constraint bounds
lc = [ lift_req, projected_area(wing)]
uc = [ lift_req, projected_area(wing)]

# Closures for in-place assignment and parameters
make_surface(wing) = make_surface(wing, span_num, chord_num)

function obj_func!(x)
    surf = (make_surface ∘ make_wing)(x[1:end-1])
    run_case!(surf, x[end], system, state)
    evaluate_induced_drag(system)
end

function cons_func!(R, x)
    wing = make_wing(x[1:end-1])
    surf = make_surface(wing)
    run_case!(surf, x[end], system, state)
    R .= evaluate_constraints(wing, system)
end

## Optim variables
optimize_chords = TwiceDifferentiable(obj_func!, x0)
cons_lift_area  = TwiceDifferentiableConstraints(cons_func!, l_bound, u_bound, lc, uc)

## Evaluate case
res_chord = optimize(optimize_chords,          # Objective functions
                     cons_lift_area,           # Constraint
                     x0,                       # Initial value
                     IPNewton(),               # Optimization algorithm
                    #  autodiff = :forward,    # Automatic differentiation
                     Optim.Options(
                                   # extended_trace = true,
                                   show_trace = true
                                  )
                    )

## Plotting
#============================================#

ini_plan = plot_wing(wing)

opt_wing = make_wing(res_chord.minimizer[1:end-1])
opt_plan = plot_wing(opt_wing)

b = span(opt_wing)
plot(xlim = (-b/2, b/2), zlim = (-b/2, b/2))
plot!(ini_plan, label = "Initial")
plot!(opt_plan, color = :blue, label = "Optimized")