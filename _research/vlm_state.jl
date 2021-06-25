## Aircraft analysis case
using Revise
using AeroMDAO
using BenchmarkTools
using ForwardDiff
using StaticArrays

## Wing
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

## Assembly
wing_panels,  wing_normals  = panel_wing(wing,  20, 10;)
htail_panels, htail_normals = panel_wing(htail, 12, 12;
                                          position = [4., 0, 0],
                                          angle    = deg2rad(-2.),
                                          axis     = [0., 1., 0.]
                                         )
vtail_panels, vtail_normals = panel_wing(vtail, 12, 10; 
                                         position = [4., 0, 0],
                                         angle    = π/2, 
                                         axis     = [1., 0., 0.]
                                        )

aircraft = Dict(
                "Wing"              => (wing_panels,  wing_normals),
                "Horizontal Tail" 	=> (htail_panels, htail_normals),
                "Vertical Tail"   	=> (vtail_panels, vtail_normals),
               );

S, b, c  = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
wing_mac = mean_aerodynamic_center(wing);

## Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [ wing_mac[1], 0., 0.]
V, α, β = 1.0, 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω);

## Pure
println("Pure -")
@time begin
    data = 
    solve_case(aircraft, fs; 
               rho_ref     = ρ, 		# Reference density
               r_ref       = ref, 		# Reference point for moments
               area_ref    = S, 		# Reference area
               span_ref    = b, 		# Reference span
               chord_ref   = c, 		# Reference chord
               name        = ac_name,	# Aircraft name
            #    print       = true,		# Prints the results for the entire aircraft
            #    print_components = true,	# Prints the results for each component
              );
end

## Data collection
comp_names = (collect ∘ keys)(data)
comp = comp_names[1]
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normies, horses, Γs = data[comp];
print_coefficients(nf_coeffs, ff_coeffs, comp)

## Impure
println("Impure -")
@time begin
    # Set up state
    state = VLMState(fs.V, fs.alpha, fs.beta, fs.omega; 
                     r_ref     = [ wing_mac[1], 0., 0 ],  
                     rho_ref   = ρ,
                     area_ref  = S, 
                     chord_ref = c, 
                     span_ref  = b,
                     name      = ac_name);

    # Solve system
    system, surfs = solve_case!(aircraft, state)
end

## Get coefficients
rate_coeffs = rate_coefficient(state)
coeffs      = aerodynamic_coefficients(surfs, state)

## Print coefficients for component of your choice
wing_names = (collect ∘ keys)(coeffs)
comp       = wing_names[1]
nf_t, ff_t = coeffs[comp]

print_coefficients(nf_t, ff_t, comp)

## Get component of your choice and its properties
surf_names = (collect ∘ keys)(surfs)
comp       = surf_names[1]
surf       = surfs[comp]

Γs     = circulations(surf)
horses = horseshoes(surf)
pts    = collocation_points(surf)
Fs     = surface_forces(surf)
Ms     = surface_moments(surf)

## Compute aerodynamic coefficients of your choice
CFs = surface_force_coefficients(surf, state)
CMs = surface_moment_coefficients(surf, state)
nf  = nearfield_coefficients(surf, state)
ff  = farfield_coefficients(surf, state)

print_coefficients(nf, ff, comp)

## VLMState Derivatives
function vlm_state_derivatives(fs)
    # Closure
    function get_derivatives(x :: AbstractVector{<: Real})
        state = VLMState(x[1], x[2], x[3], x[4:6],
                         rho_ref = ρ,
                         area_ref = S, 
                         chord_ref = c, 
                         span_ref = b,
                         name = ac_name);

        freestream_to_cartesian(-state.U, state.alpha, state.beta)
    end

    V, α, β, Ω = fs.V, fs.alpha, fs.beta, fs.omega
    jac = ForwardDiff.jacobian(get_derivatives, [ V, α, β, Ω... ])
end

dvs = vlm_state_derivatives(fs)

## VLMSystem Derivatives
function vlm_system_derivatives(system, surfs, fs)
    # Closure
    function get_derivatives(x :: AbstractVector{<: Real})
        state = VLMState(x[1], x[2], x[3], x[4:6],
                         r_ref = x[7:end],
                         rho_ref = ρ,
                         area_ref = S, 
                         chord_ref = c, 
                         span_ref = b,
                         name = ac_name);

        compute_influence_matrix!(system, state.V)

        AIC(system)
    end

    V, α, β, Ω, r_ref = fs.V, fs.alpha, fs.beta, fs.omega, [0.25, 0., 0. ]
    jac = ForwardDiff.jacobian(get_derivatives, [ V, α, β, Ω..., r_ref... ])
end

dvs = vlm_system_derivatives(system, (collect ∘ values)(surfs), fs)