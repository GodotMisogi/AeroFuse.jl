##
using AeroMDAO
using BenchmarkTools

## Aircraft definitions
#==========================================================================================#

# Wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 2.0],
            spans     = [5.0],
            dihedrals = [11.31],
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
                    sweep_LEs = [7.97])

wing_panels  = panel_wing(wing, 20, 10, spacing = "cosine")

htail_panels = panel_wing(htail, 6, 6;
                            position = [4., 0, 0],
                            angle    = deg2rad(0.),
                            axis 	   = [0., 1., 0.],
                            spacing  = "cosine"
                            )

vtail_panels = panel_wing(vtail, 6, 5; 
                            position = [4., 0, 0],
                            angle    = π/2, 
                            axis 	   = [1., 0., 0.],
                            spacing  = "cosine"
                            )

aircraft = Dict("Wing"            => wing_panels,
                "Horizontal Tail" => htail_panels,
                "Vertical Tail"   => vtail_panels)

ρ       = 1.225
x_ref   = [0.5, 0., 0.]
S, b, c = 9.0, 10.0, 0.9
αs      = -5:5

## Angle of attack sweep
#==========================================================================================#

# Stateful
function alpha_sweep!(α, system :: VLMSystem, surfs, state :: VLMState)
    state.alpha = α
    solve_case!(system, (collect ∘ values)(surfs), state)
    aerodynamic_coefficients(surfs, state)
end

# Stateful
println("AeroMDAO Aircraft Stateful -")
@time begin
    # Set up state
    state = VLMState(1., 0., 0., zeros(3); 
                     r_ref     = x_ref,  
                     rho_ref   = ρ,
                     area_ref  = S, 
                     chord_ref = c, 
                     span_ref  = b);

    # Set up system and surfaces                     
    system, surfs = build_system(aircraft);
    
    # Evaluate sweep
    coeffs = alpha_sweep!.(deg2rad.(αs), Ref(system), Ref(surfs), Ref(state))
end

## "Functional"
println("AeroMDAO Aircraft Functional -")
@time begin
    # Build vector of Freestreams
    fses = Freestream.(1.0, αs, 0., Ref(zeros(3)))

    # Evaluate sweep
    data = solve_case.(Ref(aircraft), fses; 
                       rho_ref   = ρ,
                       r_ref     = x_ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c)
end
