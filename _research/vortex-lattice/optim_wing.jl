## Optimisation tests
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
    CDi  = ff[1]

    lift, area, CDi
end

## Optimization
#============================================#

# Parameters
V, α, ρ  = 29., 5., 1.225

# Design variables
n    = 4
wing = Wing(foils     = fill(Foil(naca4(2,4,1,2)), n),
            chords    = fill(0.314, n),
            twists    = fill(0.0, n),
            spans     = fill(1.3/(n-1), n-1),
            dihedrals = fill(0., n-1),
            LE_sweeps = fill(0., n-1))

x0 = [(chords ∘ right)(wing); α]

make_wing(x) = Wing(chords    = x,
                    foils     = foils(right(wing)),
                    spans     = spans(right(wing)),
                    twists    = rad2deg.(twists(right(wing))),
                    dihedrals = rad2deg.(dihedrals(right(wing))),
                    LE_sweeps = rad2deg.(sweeps(right(wing))))


# wing = WingSection(root_foil  = naca4((2,4,1,2)),
#                    tip_foil   = naca4((2,4,1,2)),
#                    span       = 2.6,
#                    root_chord = 0.314,
#                    taper      = 0.8,
#                    dihedral   = 5.0,
#                    LE_sweep   = 10.)

# x0 = [ mean_aerodynamic_chord(wing); taper_ratio(right(wing)); α ]

# make_wing(x) = WingSection(root_foil  = naca4((2,4,1,2)),
#                            tip_foil   = naca4((2,4,1,2)),
#                            span       = span(wing),
#                            root_chord = x[1],
#                            taper      = x[2],
#                            dihedral   = rad2deg((dihedrals ∘ right)(wing)[1]),
#                            LE_sweep   = rad2deg((sweeps ∘ right)(wing)[1]))


# Meshing and assembly
wing_mac = mean_aerodynamic_center(wing);
b, S, c  = info(wing)[1:end-1]

span_num  = 10
chord_num = 1

# Test
test_CDi = @time evaluate_CDi(wing, V, α, ρ, span_num, chord_num)
test_con = @time evaluate_cons(wing, V, α, ρ, span_num, chord_num)

# Bounds and constraints
weight      = 12 * 9.81
load_factor = 1.5
lift_req    = weight * load_factor

l_bound = fill(1e-12, length(x0))
u_bound = fill(Inf,   length(x0))

lc = [ lift_req, projected_area(wing), 0. ]
uc = [ lift_req, projected_area(wing), Inf ]

# Closures
evaluate_CDi(x) = evaluate_CDi(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)
cons!(cs, x) = cs .= evaluate_cons(make_wing(x[1:end-1]), V, x[end], ρ, span_num, chord_num)

# Optim variables
optimize_chords = TwiceDifferentiable(evaluate_CDi, x0)
cons_lift_area  = TwiceDifferentiableConstraints(cons!, l_bound, u_bound, lc, uc);

## Run optimisation
res_func = optimize(optimize_chords,      # Objective functions
                    cons_lift_area,       # Constraint
                    x0,                   # Initial value
                    IPNewton(),           # Optimization algorithm
                    autodiff = :forward,  # Automatic differentiation
                    Optim.Options(
                                #   extended_trace = true,
                                  show_trace = true
                                 )
                   )

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

## STATEFUL TEST
#==========================================================================================#

# Helper functions
#============================================#

function make_surface(wing, name, span_num, chord_num)
    wing_panels, wing_normals = panel_wing(wing, span_num, chord_num;)
    VLMSurface(wing_panels, wing_normals, name)
end

function run_case!(surf, α, system, state)
    # Update angle of attack and surface
    state.alpha        = α
    system.surfaces[findfirst(x -> name(surf) == x, name.(surfaces(system)))] = surf

    # Run case
    evaluate_case!(system, state)
end

# Objective function
evaluate_induced_drag(system, state) = sum(x -> farfield_coefficients(x, state), surfaces(system))[1]

# Lift and area constraint function
evaluate_constraints(wing, system) = [ sum(farfield_forces, surfaces(system))[3]; projected_area(wing); sum(x -> farfield_coefficients(x, state), surfaces(system))[1] ]

## Optimisation
#============================================#

# Set up state
state = VLMState(V, deg2rad(α), 0., [0., 0., 0.],
                 rho_ref   = ρ,
                 r_ref     = [ wing_mac[1], 0, 0 ],
                 area_ref  = projected_area(wing),
                 chord_ref = mean_aerodynamic_chord(wing),
                 span_ref  = span(wing));

# Aircraft assembly
htail   = WingSection(root_foil  = naca4((0,0,1,2)),
                      tip_foil   = naca4((0,0,1,2)),
                      span       = 0.87,
                      taper      = 1.0,
                      root_chord = 0.215)

vtail_1 = HalfWingSection(root_foil  = naca4((0,0,1,2)),
                          tip_foil   = naca4((0,0,1,2)),
                          span       = 0.4,
                          taper      = 0.731,
                          root_chord = 0.260,
                          LE_sweep   = 15.)
vtail_2  = vtail_1
aircraft = Dict("Wing"    => panel_wing(wing, span_num, chord_num),
                "HTail"   => panel_wing(htail, 6, 6;
                                        position = [1., 0, 0.],
                                        angle    = deg2rad(0.),
                                        axis     = [0., 1., 0.]),
                # "VTail 1" =>  panel_wing(vtail_1, 6, 5;
                #                          position = [1., 0., 0],
                #                          angle    = π/2,
                #                          axis     = [1., 0, 0],
                #                          spacing  = "cosine"),
                # "VTail 2" => panel_wing(vtail_2, 6, 5;
                #                         position = [1., 0.435, 0],
                #                         angle    = π/2,
                #                         axis     = [1., 0, 0],
                #                         spacing  = "cosine")
                );

# Solve initial case
system   = solve_case(aircraft, state);

# Test case - Fixed speed
state.speed   = V
state.alpha   = deg2rad(α)
state.rho_ref = ρ

# Closures for in-place assignment and parameters
make_surface(wing) = make_surface(wing, "Wing", span_num, chord_num)

function obj_func!(x)
    wing = make_wing(x[1:end-1])
    surf = make_surface(wing)
    state.area_ref = projected_area(wing)
    run_case!(surf, x[end], system, state)
    evaluate_induced_drag(system, state)
end

function cons_func!(R, x)
    wing = make_wing(x[1:end-1])
    surf = make_surface(wing)
    state.area_ref = projected_area(wing)
    run_case!(surf, x[end], system, state)
    R .= evaluate_constraints(wing, system)
end

## Optim variables
optimize_chords = TwiceDifferentiable(obj_func!, x0)
cons_lift_area  = TwiceDifferentiableConstraints(cons_func!, l_bound, u_bound, lc, uc);

## Evaluate case
res_state = optimize(optimize_chords,          # Objective functions
                     cons_lift_area,           # Constraint
                     x0,                       # Initial value
                     IPNewton(),               # Optimization algorithm
                    #  autodiff = :forward,    # Automatic differentiation
                     Optim.Options(
                                   # extended_trace = true,
                                   show_trace = true
                                  )
                    )

## Sanity check
x_opt    = res_state.minimizer
CDi_opt  = evaluate_CDi(x_opt)
cons_opt = cons!(zeros(3), x_opt)

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