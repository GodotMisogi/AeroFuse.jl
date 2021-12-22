## Aircraft analysis case
using Revise
using AeroMDAO
using ForwardDiff
using LinearAlgebra

## Lifting surfaces setup
#==========================================================================================#

# Wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 2.0],
            spans     = [5.0],
            dihedrals = [11.3],
            LE_sweeps = [2.29]);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97]);

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


## Case
ac_name  = "My Aircraft"
S, b, c  = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
wing_mac = mean_aerodynamic_center(wing);
ρ        = 1.225
x_w      = [ wing_mac[1], 0., 0 ]
fs       = Freestream(1.0, 1.0, 0.0, [0.0, 0.0, 0.0])

@time begin
    # Set up state
    state = VLMState(fs;
                     r_ref     = x_w,
                     rho_ref   = ρ,
                     area_ref  = S,
                     chord_ref = c,
                     span_ref  = b,
                     name      = ac_name);

    # Solve system
    system = solve_case(aircraft, state)
end;

## Get surfaces and coefficients
surfs       = surfaces(system)
rate_coeffs = rate_coefficient(state)
coeffs      = aerodynamic_coefficients(surfs, state)

## Print coefficients for all surfaces
print_coefficients(surfs, state);

## Print coefficients for component of your choice
print_coefficients(surfs[1], state)

## Get component of your choice and its properties
surf   = surfs[1]

println(name(surf))

Γs     = circulations(surf)
horses = horseshoes(surf)
pts    = collocation_points(surf)
Fs     = surface_forces(surf)
Ms     = surface_moments(surf)

## Compute aerodynamic coefficients of your choice
CFs = surface_force_coefficients(surf, state)
CMs = surface_moment_coefficients(surf, state)
nf  = nearfield_coefficients(surf, state)
ff  = farfield_coefficients(surf, state);

print_coefficients(nf, ff, name(surf))

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

        freestream_to_cartesian(-state.speed, state.alpha, state.beta)
    end

    V, α, β, Ω = fs.speed, fs.alpha, fs.beta, fs.omega
    jac = ForwardDiff.jacobian(get_derivatives, [ V; α; β; Ω ])
end

dvs = vlm_state_derivatives(fs)

## VLMSystem Derivatives
function vlm_system_derivatives(system, surf, fs, r_ref)
    # Closure
    function get_derivatives(x)
        state = VLMState(x[1], x[2], x[3], x[4:6],
                         r_ref = x[7:end],
                         rho_ref = ρ,
                         area_ref = S,
                         chord_ref = c,
                         span_ref = b,
                         name = ac_name);

        Γs = AeroMDAO.VortexLattice.solve_system(horseshoes(system), normals(system), normalize(state.velocity), state.omega)

        # system.circulations = Γs
        # update_circulations!(system)

        forces, moments, trefftz_force = AeroMDAO.VortexLattice.evaluate_dynamics(Γs, horseshoes(system), state.velocity, state.alpha, state.beta, state.omega, state.rho_ref, state.r_ref)

        # nearfield_coeffs, farfield_coeffs, CFs, CMs = AeroMDAO.VortexLattice.evaluate_coefficients(forces, moments, trefftz_force, state.velocity, state.alpha, state.beta, state.rho_ref, state.area_ref, state.chord_ref, state.span_ref)

        [ nearfield_coeffs[:]; farfield_coeffs[:]; reduce(vcat, CFs[:]); reduce(vcat, CMs[:]) ]
    end

    V, α, β, Ω = fs.speed, fs.alpha, fs.beta, fs.omega
    jac = ForwardDiff.jacobian(get_derivatives, [ V; α; β; Ω; r_ref ])
end

dvs = vlm_system_derivatives(system, surf, fs, x_w)