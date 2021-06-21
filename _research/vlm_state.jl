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
wing_panels,  wing_normals  = panel_wing(wing,  20, 10);
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
                "Wing" 			  	=> (wing_panels,  wing_normals),
                "Horizontal Tail" 	=> (htail_panels, htail_normals),
                "Vertical Tail"   	=> (vtail_panels, vtail_normals),
               );

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
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
end;

## Data collection
comp_names = (collect ∘ keys)(data)
comp = comp_names[1]
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horses, Γs = data[comp];
print_coefficients(nf_coeffs, ff_coeffs, comp)

## Impure
println("Impure -")
@time begin
    # Set up state
    state = VLMState(fs, 
                     rho_ref   = ρ,
                     r_ref     = SVector(ref...),
                     area_ref  = S, 
                     chord_ref = c, 
                     span_ref  = b);

    # Solve system
    system, surfs, coeffs = solve_case!(aircraft, state);
end;

## Get only coefficients
wing_names = (collect ∘ keys)(coeffs)
comp       = wing_names[2]
nf_t, ff_t = coeffs[comp]

print_coefficients(nf_t, ff_t, comp)

## Get variables of your choice
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

## Derivatives
function wtf(aircraft, fs)
    
    # Closure
    function get_derivatives(x)
        state = VLMState(fs.V, x[1], x[2], x[3:end], 
                         rho_ref   = ρ,
                         r_ref     = ref,
                         area_ref  = S, 
                         chord_ref = c, 
                         span_ref  = b, 
                         name      = ac_name);

        system, surfs, nf_t, ff_t = AeroMDAO.VortexLattice.solve_case!(aircraft, state);

        nf_t
    end

    jac = ForwardDiff.jacobian(get_derivatives, [fs.alpha, fs.beta, fs.omega... ])
end

dvs = wtf(aircraft, fs)

##
