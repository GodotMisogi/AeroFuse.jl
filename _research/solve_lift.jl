## Nonlinear solution cases
using Revise
using AeroMDAO
using NLsolve
using StaticArrays
using BenchmarkTools

## System setup
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = WingSection(root_foil  = naca4((4,4,1,2)),
                   span       = 1.3,
                   dihedral   = 0.0,
                   sweep_LE   = 15.0,
                   taper      = 1.0,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
name        = "Wing"
print_info(wing)

# Mesh
span_num        = 20
chord_num       = 13
panels, normies = panel_wing(wing, span_num, chord_num);
aircraft        = Dict(name => (panels, normies));

## Define VLMState and VLMSystem

# Set up state
V, α, β = 20., 0., 0.
state = VLMState(Freestream(V, α, β, [0.0, 0.0, 0.0]), 
                 rho_ref   = 1.225,
                 r_ref     = SVector(wing_mac[1], 0., 0.),
                 area_ref  = projected_area(wing), 
                 chord_ref = mean_aerodynamic_chord(wing), 
                 span_ref  = span(wing));

# Build system
system, surfs = build_system(aircraft);

## Weight variables
weight      = 15 * 9.81
load_factor = 1.5;

## Residual setup
#==========================================================================================#

load_factor_residual(L, W, n) = L - n * W

function freestream_case!(system :: VLMSystem, state :: VLMState, surfs, weight, load_factor)
    V = freestream_to_cartesian(-state.U, state.alpha, state.beta)

    generate_system!(system, V, state.omega) # Pre-allocated version for efficiency
    solve_system!(system)
    update_circulations!(circulations(system), surfs)

    # Compute lift
    compute_surface_forces!.(surfs, Ref(system), Ref(V), Ref(state.omega), state.rho_ref)
    lift = body_to_wind_axes(sum(sum ∘ surface_forces, surfs), state.alpha, state.beta)[3]

    # Weight residual
    load_factor_residual(lift, weight, load_factor)
end

## Angle of attack residual
#==========================================================================================#

function solve_alpha_residual!(R, x, system :: VLMSystem, state :: VLMState, surfs, weight, load_factor)
    state.alpha = x[1]
    R .= freestream_case!(system, state, surfs, weight, load_factor)
end

## Test case - Fixed speed
state.U = V
x       = [ state.alpha ]

# Nonlinear solution
solve_alpha_residual!(R, x) = solve_alpha_residual!(R, x, system, state, values(surfs), weight, load_factor)
@time res_alpha = 
    nlsolve(solve_alpha_residual!, x,
            method   = :newton,
            # autodiff = :forward, # NOT WORKING
            # show_trace = true
           )

## Check numbers
state.alpha          = res_alpha.zero[1]
system, surfs        = solve_case!(aircraft, state)
coeffs               = aerodynamic_coefficients(surfs, state)
L                    = force(coeffs["Wing"][2][1:3], dynamic_pressure(state.rho_ref, state.U), state.area_ref)[3]
n_calc               = L / weight

println("Load factor: $n_calc")
println("Weight: $weight, N")
println("Lift: $L, N")
println("Speed: $(state.U), m/s")
println("Angle of attack: $(rad2deg(res_alpha.zero[1]))ᵒ")

## Speed residual
#==========================================================================================#

function solve_speed_residual!(R, x, system :: VLMSystem, state :: VLMState, surfs, weight, load_factor)
    state.U = x[1]
    R .= freestream_case!(system, state, surfs, weight, load_factor)
end

# Test case - Fixed angle of attack
state.alpha = deg2rad(3.)
x           = [ state.U ]

# Nonlinear solution
solve_speed_residual!(R, x) = solve_speed_residual!(R, x, system, state, values(surfs), weight, load_factor)
@time res_speed = 
    nlsolve(solve_speed_residual!, x,
            method   = :newton,
            # autodiff = :forward, # NOT WORKING
            # show_trace = true,
           )

## Check numbers
state.U              = res_speed.zero[1]
system, surfs        = solve_case!(aircraft, state)
coeffs               = aerodynamic_coefficients(surfs, state)
L                    = force(coeffs["Wing"][2][1:3], dynamic_pressure(state.rho_ref, state.U), state.area_ref)[3]
n_calc               = L / weight

println("Load factor: $n_calc")
println("Weight: $weight, N")
println("Lift: $L, N")
println("Speed: $(res_speed.zero[1]), m/s")
println("Angle of attack: $(rad2deg(state.alpha))ᵒ")