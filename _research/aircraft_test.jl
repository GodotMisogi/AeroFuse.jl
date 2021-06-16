##
using AeroMDAO

## Wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [0.0, 0.0],
            spans     = [5.0],
            dihedrals = [11.39],
            sweep_LEs = [0.]);
print_info(wing, "Wing")

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
htail = Wing(foils     = htail_foils,
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39])
print_info(htail, "Horizontal Tail")

# Vertical tail
vtail_foils = Foil.(fill(naca4((0,0,0,9)), 2))
vtail = HalfWing(foils     = vtail_foils, 
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97])
print_info(vtail, "Vertical Tail")

## Assembly
wing_panels  = panel_wing(wing, [25], 15;
                          spacing = "cosine"
                         )
htail_panels = panel_wing(htail, [6], 6;
                          position	= [4., 0, 0],
                          angle 	= deg2rad(-2.),
                          axis 	  	= [0., 1., 0.],
                          spacing   = "uniform"
                         )
vtail_panels = panel_wing(vtail, [5], 6; 
                          position 	= [4., 0, 0],
                          angle 	= π/2, 
                          axis 	 	= [1., 0., 0.],
                          spacing   = "uniform"
                         )

aircraft = Dict("Wing" 			  	=> wing_panels,
                "Horizontal Tail" 	=> htail_panels,
                "Vertical Tail"   	=> vtail_panels)

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

## Stability case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, α, β = 1.0, 1.0, 1.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

dv_data = 
    solve_stability_case(aircraft, fs;
                         rho_ref     = ρ,
                         r_ref       = ref,
                         area_ref    = S,
                         span_ref    = b,
                         chord_ref   = c,
                         name        = ac_name,
                         print       = true,
                         print_components = true,
                        );

## Data collection
names = (collect ∘ keys)(dv_data) 
comp  = names[1]
nf, ff, dvs = dv_data[comp];
print_coefficients(nf, ff, comp)
print_derivatives(dvs, comp)
