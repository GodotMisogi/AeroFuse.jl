## Analysis for Redbing Bloodstream
using Revise
using AeroMDAO
using NLsolve

## System setup
#==========================================================================================#

# Define lifting surfaces
wing = WingSection(root_foil  = naca4((2,4,1,2)),
                   tip_foil   = naca4((2,4,1,2)),
                   span       = 1.3,
                   taper      = 1.0,
                   root_chord = 0.314)
wing_name = "Wing"
wing_mac = mean_aerodynamic_center(wing)

htail = WingSection(root_foil  = naca4((0,0,1,2)),
                    tip_foil   = naca4((0,0,1,2)),
                    span       = 0.435,
                    taper      = 1.0,
                    root_chord = 0.215)
htail_name = "Horizontal Tail"

vtail_1 = HalfWingSection(root_foil  = naca4((0,0,1,2)),
                          tip_foil   = naca4((0,0,1,2)),
                          span       = 0.15,
                          taper      = 0.731,
                          root_chord = 0.260)
vtail_1_name = "Vertical Tail 1"

vtail_2 = vtail_1 
vtail_2_name = "Vertical Tail 2"


print_info(wing, wing_name)
print_info(htail, htail_name)
print_info(vtail_1, vtail_1_name)
print_info(vtail_2, vtail_2_name)

# Mesh
wing_panels, wing_normals       = panel_wing(wing,  20, 13;
                                             angle = deg2rad(1.9),
                                             axis  = [0., 1., 0.])
htail_panels, htail_normals     = panel_wing(htail, 12, 12;
                                             position = [1., 0, 0.15],
                                             angle    = deg2rad(0.),
                                             axis     = [0., 1., 0.])
vtail_1_panels, vtail_1_normals = panel_wing(vtail_1, 6, 5; 
                                             position = [1., -0.435, 0],
                                             angle    = π/2, 
                                             axis     = [1., 0, 0],
                                             spacing  = "cosine")
vtail_2_panels, vtail_2_normals = panel_wing(vtail_2, 6, 5; 
                                             position = [1., 0.435, 0],
                                             angle    = π/2, 
                                             axis     = [1., 0, 0],
                                             spacing  = "cosine")

aircraft = Dict(
                wing_name    => (wing_panels   , wing_normals   ),
                htail_name   => (htail_panels  , htail_normals  ),
                vtail_1_name => (vtail_1_panels, vtail_1_normals),
                vtail_2_name => (vtail_2_panels, vtail_2_normals)
                );

## System-state evaluation
state = VLMState(Freestream(20., 1., 0., [0., 0., 0.]), 
                 rho_ref   = 1.225,
                 r_ref     = [ wing_mac[1], 0, 0 ],
                 area_ref  = projected_area(wing), 
                 chord_ref = mean_aerodynamic_chord(wing), 
                 span_ref  = span(wing))

system, surfs = solve_case(aircraft, state);

## Printing coefficients
print_coefficients(surfs, state);

## Residual setup
#==========================================================================================#

load_factor_residual(L, W, n) = L - n * W

function load_factor_case!(system :: VLMSystem, surfs :: Vector{<: VLMSurface}, state :: VLMState, weight, load_factor)
    # Evaluate aerodynamics
    evaluate_case!(system, surfs, state)

    # Compute lift
    lift = body_to_wind_axes(sum(sum ∘ surface_forces, surfs), state.alpha, state.beta)[3]

    # Weight residual
    load_factor_residual(lift, weight, load_factor)
end

## Angle of attack residual
#==========================================================================================#

function solve_alpha_residual!(R, x, system :: VLMSystem, surfs :: Vector{<: VLMSurface}, state :: VLMState, weight, load_factor)
    state.alpha = x[1]
    R .= load_factor_case!(system, surfs, state, weight, load_factor)
end

## Test case - Fixed speed
system, surfs = build_system(aircraft)
weight        = 15 * 9.81
load_factor   = 1.0
state.speed       = 25.
state.rho_ref = 0.98
x             = [ state.alpha ]

# Nonlinear solution
solve_alpha_residual!(R, x) = solve_alpha_residual!(R, x, system, surfs, state, weight, load_factor)
@time res_alpha = 
    nlsolve(solve_alpha_residual!, x,
            method     = :newton,
            show_trace = true,
            # autodiff = :forward, # NOT WORKING
            # extended_trace = true
           )

## Check numbers
lift     = sum(sum ∘ surface_forces, values(surfs))[3]
load_fac = lift / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(state.speed) m/s")
println("Angle of attack: $(rad2deg(state.alpha))ᵒ")

## Speed residual
#==========================================================================================#

function solve_speed_residual!(R, x, system :: VLMSystem, surfs :: Vector{<: VLMSurface}, state :: VLMState, weight, load_factor)
    state.speed = x[1]
    R .= load_factor_case!(system, surfs, state, weight, load_factor)
end

# Test case - Fixed angle of attack
weight        = 15 * 9.81
load_factor   = 1.5
state.alpha   = deg2rad(3.)
state.rho_ref = 1.1
x             = [ state.speed ]

# Nonlinear solution
solve_speed_residual!(R, x) = solve_speed_residual!(R, x, system, surfs, state, weight, load_factor)
@time res_speed = 
    nlsolve(solve_speed_residual!, x,
            method     = :newton,
            show_trace = true,
            # autodiff = :forward, # NOT WORKING
            # extended_trace = true
           )

## Check numbers
lift     = sum(sum ∘ surface_forces, values(surfs))[3]
load_fac = lift / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(state.speed) m/s")
println("Angle of attack: $(rad2deg(state.alpha))ᵒ")

## Plotting
#==========================================================================================#

function plot_aircraft(aircraft)
    panels = (vec ∘ first).(values(aircraft))
    reduce(vcat, plot_panels.(panels))
end

pans = plot_aircraft(aircraft)

using Plots

plot(aspect_ratio = 1, zlim = (-span(wing) / 2, span(wing) / 2))
plot!.(pans, color = :black, label = :none)
plot!()