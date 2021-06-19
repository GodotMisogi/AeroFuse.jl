##
using AeroMDAO
using NLsolve

## Stability analysis
function aircraft_nearfield_forces(aircraft, fs, ρ, ref, S, b, c)
    nf, ff, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = 
        solve_case(aircraft, fs; 
                   rho_ref     = ρ, 
                   r_ref       = ref, 
                   area_ref    = S, 
                   span_ref    = b, 
                   chord_ref   = c,
                   )["Aircraft"]
    
    nf
end

pitching_moment_coefficient(aircraft, α :: Real, ρ, ref, S, b, c) = aircraft_nearfield_forces(aircraft, Freestream(1.0, α, 0., zeros(3)), ρ, ref, S, b, c)[5]

## Define aircraft

# Wing
wing  = WingSection(span       = 4.0,
                    dihedral   = 5.0,
                    sweep_LE   = 15.0,
                    taper      = 0.4,
                    root_chord = 2.0,
                    root_twist = 0.0,
                    tip_twist  = -2.0,
                    root_foil  = naca4((2,4,1,2)),
                    tip_foil   = naca4((2,4,1,2)))
wing_mac  = mean_aerodynamic_center(wing)
wing_pos  = [0., 0., 0.]
x_w       = [ wing_pos[1]; 0.; 0. ]
wing_plan = plot_wing(wing;  
                      position = wing_pos)

print_info(wing, "Wing")

htail = WingSection(span       = 1.0,
                    dihedral   = 0.0,
                    sweep_LE   = 15.0,
                    taper      = 0.6,
                    root_chord = 0.8,
                    root_twist = 0.0,
                    tip_twist  = 0.0,
                    root_foil  = naca4((0,0,1,2)),
                    tip_foil   = naca4((0,0,0,9)));
htail_mac  = mean_aerodynamic_center(htail)
htail_pos  = [5., 0., 0.]
α_h_i      = 0.
htail_plan = plot_wing(htail;
                       position = htail_pos)

print_info(htail, "Horizontal Tail")

vtail = HalfWingSection(span       = 0.8,
                        dihedral   = 0.0,
                        sweep_LE   = 8.0,
                        taper      = 0.6,
                        root_chord = 0.8,
                        root_twist = 0.0,
                        tip_twist  = 0.,
                        root_foil  = naca4((0,0,0,9)),
                        tip_foil   = naca4((0,0,0,9)))
vtail_mac  = mean_aerodynamic_center(vtail) # NEEDS FIXING FOR ROTATION
vtail_pos  = [5., 0., 0.]
vtail_plan = plot_wing(vtail;
                       position = vtail_pos,
                       angle    = π/2)
                        
## Panelling and assembly
wing_panels  = panel_wing(wing, [20], 10;
                          position = wing_pos
                         )
htail_panels = panel_wing(htail, [10], 5;
                          position = htail_pos,
                          angle    = deg2rad(α_h_i),
                          axis 	   = [0., 1., 0.]
                         )
vtail_panels = panel_wing(vtail, [10], 5;
                          position = vtail_pos,
                          angle    = π/2
                         )

aircraft     = Dict("Wing" 			  => wing_panels,
                    "Horizontal Tail" => htail_panels,
                    "Vertical Tail"   => vtail_panels);

## Evaluate case
ac_name = "My Aircraft"
ρ 		= 1.225
ref     = x_w
V, α, β = 1.0, 0.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 	    = Freestream(V, α, β, Ω)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

## Evaluate trim
trim_state!(R, α) = R .= pitching_moment_coefficient(aircraft, α[1], ρ, ref, S, b, c)
prob = nlsolve(trim_state!, [-1.0], 
               method = :newton,
               autodiff = :forward,
               show_trace = true, 
            #    extended_trace = true
               )
