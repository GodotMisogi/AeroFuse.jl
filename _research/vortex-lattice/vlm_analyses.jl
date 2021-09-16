## Nonlinear solution cases
using Revise
using AeroMDAO
using NLsolve
using DataFrames

## Lifting surfaces setup
#==========================================================================================#

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 2.0],
            spans     = [5.0],
            dihedrals = [11.3],
            sweep_LEs = [2.29]);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97]);

## Meshing and assembly
wing_panels,  wing_normals  = panel_wing(wing,  20, 13;)
htail_panels, htail_normals = panel_wing(htail, 6, 6;
                                         position = [4., 0, 0],
                                         angle    = deg2rad(-2.),
                                         axis     = [0., 1., 0.]
                                        )
vtail_panels, vtail_normals = panel_wing(vtail, 6, 5;
                                         position = [4., 0, 0],
                                         angle    = π/2,
                                         axis     = [1., 0., 0.]
                                        )

# Aircraft assembly
aircraft = Dict(
                "Wing"            => (wing_panels,  wing_normals),
                "Horizontal Tail" => (htail_panels, htail_normals),
                "Vertical Tail"   => (vtail_panels, vtail_normals),
               );

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

## Angle of attack sweep
#==========================================================================================#

function alpha_sweep!(α, system, state)
    state.alpha = α
    evaluate_case!(system, state)
    nf, ff = aerodynamic_coefficients(surfaces(system), state)[state.name]

    [ ff; nf[4:end] ]
end

##
αs = deg2rad.(-5:5)
@time results = alpha_sweep!.(αs, Ref(system), Ref(state))
#
data = DataFrame([ (α, CD, CY, CL, Cl, Cm, Cn) for (α, (CD, CY, CL, Cl, Cm, Cn)) in zip(rad2deg.(αs), results) ][:])
rename!(data, [:α, :CD, :CY, :CL, :Cl, :Cm, :Cn])

## Residual setup
#==========================================================================================#

load_factor_residual(L, W, n) = L - n * W

function load_factor_case!(system, state, weight, load_factor)
    # Evaluate aerodynamics
    evaluate_case!(system, state)

    # Compute lift
    lift = body_to_wind_axes(sum(sum ∘ surface_forces, surfaces(system)), state.alpha, state.beta)[3]

    # Weight residual
    load_factor_residual(lift, weight, load_factor)
end

## Angle of attack residual
#==========================================================================================#

function solve_alpha_residual!(R, x, system, state, weight, load_factor)
    state.alpha = x[1]
    R .= load_factor_case!(system, state, weight, load_factor)
end

## Test case - Fixed speed
weight        = 15 * 9.81
load_factor   = 1.2
state.speed   = 25.
state.rho_ref = 0.98
x             = [ state.alpha ]

# Nonlinear solution
solve_alpha_residual!(R, x) = solve_alpha_residual!(R, x, system, state, weight, load_factor)
@time res_alpha =
    nlsolve(solve_alpha_residual!, x,
            method     = :newton,
            show_trace = true,
            # autodiff = :forward, # NOT WORKING
            # extended_trace = true
           )

## Check numbers
lift     = sum(sum ∘ surface_forces, surfaces(system))[3]
load_fac = lift / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(state.speed) m/s")
println("Angle of attack: $(rad2deg(state.alpha))ᵒ")

## Speed residual
#==========================================================================================#

function solve_speed_residual!(R, x, system, state, weight, load_factor)
    state.speed = x[1]
    R .= load_factor_case!(system, state, weight, load_factor)
end

# Test case - Fixed angle of attack
weight        = 12 * 9.81
load_factor   = 1.5
state.alpha   = deg2rad(3.)
state.rho_ref = 1.1
x             = [ state.speed ]

# Nonlinear solution
solve_speed_residual!(R, x) = solve_speed_residual!(R, x, system, state, weight, load_factor)
@time res_speed =
    nlsolve(solve_speed_residual!, x,
            method     = :newton,
            show_trace = true,
            # autodiff = :forward, # NOT WORKING
            # extended_trace = true
           )

## Check numbers
lift     = sum(sum ∘ surface_forces, surfaces(system))[3]
load_fac = lift / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(state.speed) m/s")
println("Angle of attack: $(rad2deg(state.alpha))ᵒ")