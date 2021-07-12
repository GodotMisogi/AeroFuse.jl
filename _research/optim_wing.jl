## 
using Revise
using StaticArrays
using AeroMDAO
using Optim

# Helper functions
#============================================#

make_wing(chords, twists, spans, dihedrals, sweeps) = Wing(chords = chords, twists = twists, spans = spans, dihedrals = dihedrals, sweep_LEs = sweeps)

function run_case(wing, V, α)
    uniform = Freestream(V, α, 0.0, zeros(3))
    ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
    nf_coeffs, ff_coeffs, horseshoe_panels, normals, horseshoes, Γs = solve_case(wing, uniform, rho_ref = 1.225, r_ref = ref, span_num = 20, chord_num = 13)

    ff_coeffs
end

# Objective function
evaluate_CDi(wing, V, α) = run_case(wing, V, α)[1]

# Lift and area constraint function
function cons_prob(wing, V, α, ρ = 1.225)
    S    = projected_area(wing)
    L    = dynamic_pressure(ρ, V) * S * run_case(wing, V, α)[3] 

    L, S
end

## Test runs
#============================================#

# Parameters
wing = Wing(chords    = [1.0, 0.5, 0.4, 0.3, 0.2],
            twists    = [0., 0., 0., 0., 0.],
            spans     = fill(2.6 / 4, 4),
            dihedrals = [0., 0., 0., 0.],
            sweep_LEs = [0., 0., 0., 0.])

V         = 20.
α         = 3.

# Test case
ff_coeffs = @time run_case(wing, V, α)
cons_test = @time cons_prob(wing, V, α, 1.225)

## Optimisation
#============================================#

# Design variables - chords and angle of attack
x0           = [chords(wing); α]
weight       = 12 * 9.81
load_factor  = 1.5
lift_req     = weight * load_factor

# Bounds
l_bound      = fill(1e-12, length(x0))
u_bound      = fill(Inf, length(x0))
bound_chords = TwiceDifferentiableConstraints(l_bound, u_bound)

# Constraint bounds
lc           = [ lift_req, projected_area(wing) ]
uc           = [ lift_req, projected_area(wing) ]

# Closures
make_wing(x)      = make_wing(x, twists(wing), spans(wing), dihedrals(wing), sweeps(wing))
evaluate_CDi(x)   = evaluate_CDi(make_wing(x[1:end-1]), V, x[end])
cons_prob!(cs, x) = cs .= cons_prob(make_wing(x[1:end-1]), V, x[end])

# Optim variables
optimize_chords = TwiceDifferentiable(evaluate_CDi, x0)
lift_constraint = TwiceDifferentiableConstraints(cons_prob!, l_bound, u_bound, lc, uc)

# Run optimisation
res_chord = optimize(optimize_chords,          # Objective functions
                     lift_constraint,          # Constraint
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

using Plots
pyplot(dpi = 300)

opt_wing  = make_wing(res_chord.minimizer[1:end-1], twists(wing), spans(wing), dihedrals(wing), sweeps(wing))
wing_plan = plot_wing(opt_wing)

b = span(opt_wing)
plot(xlim = (-b/2, b/2), zlim = (-b/2, b/2))
plot!(wing_plan, color = :blue, label = :none)


## STATEFUL TEST
#==========================================================================================#

# Lifting surfaces setup
#==========================================================================================#

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 2.0],
            spans     = [5.0],
            dihedrals = [11.3],
            sweep_LEs = [2.29]);

## Meshing and assembly
wing_panels, wing_normals = panel_wing(wing,  20, 13;)
aircraft = Dict("Wing" => (wing_panels,  wing_normals));

## System setup
#==========================================================================================#

# Set up state
wing_mac = mean_aerodynamic_center(wing);
state = VLMState(1., 0., 0., [0., 0., 0.], 
                 rho_ref   = 1.225,
                 r_ref     = [ wing_mac[1], 0, 0 ],
                 area_ref  = projected_area(wing), 
                 chord_ref = mean_aerodynamic_chord(wing), 
                 span_ref  = span(wing));

# Set up system
system = solve_case(aircraft, state);

# Helper functions
#============================================#

function make_surface(bing, span_num, chord_num)
    wing_panels, wing_normals = panel_wing(bing, span_num, chord_num;)
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

# Test case - Fixed speed
weight        = 12 * 9.81
load_factor   = 1.5
state.alpha   = deg2rad(3.)
state.rho_ref = 1.1

# Design variables - chords and angle of attack
c0        = chords(wing)
α0        = state.alpha
x0        = [c0; α0]
span_num  = 20
chord_num = 13

# Bounds
l_bound = fill(1e-12, length(x0))
u_bound = fill(Inf, length(x0))

# Constraint bounds
lc = [ load_factor * weight, projected_area(wing)]
uc = [ load_factor * weight, projected_area(wing)]

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
# bound_chords    = TwiceDifferentiableConstraints(l_bound, u_bound)

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
