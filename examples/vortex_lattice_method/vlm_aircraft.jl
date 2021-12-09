## Aircraft analysis case
using AeroMDAO
using ComponentArrays
using TimerOutputs
## Lifting surfaces

# Wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [2.0, 0.0, -2.0],
            spans     = [4.0, 0.2],
            dihedrals = [5., 30.],
            LE_sweeps = [5., 30.]);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## Assembly
wing_panels, wing_normals  = panel_wing(wing,                 # Wing or HalfWing type
                                        [20, 3],              # Number of spanwise panels for half-wing
                                        10;                   # Chordwise panels
                                        # spacing = "cosine"  # Spacing distribution: Default works well for symmetric
                                       )

htail_panels, htail_normals = panel_wing(htail, [6], 6;
                                         spacing  = Uniform()
                                        )

vtail_panels, vtail_normals = panel_wing(vtail, [6], 5;
                                         spacing  = Uniform()
                                        )

wing_horses  = Horseshoe.(wing_panels,  wing_normals)
htail_horses = Horseshoe.(htail_panels, htail_normals)
vtail_horses = Horseshoe.(vtail_panels, vtail_normals)

aircraft_panels = ComponentArray(
                                 wing  = wing_panels,
                                 htail = htail_panels, 
                                 vtail = vtail_panels
                                )

aircraft = ComponentArray(
                          wing  = wing_horses,
                          htail = htail_horses,
                          vtail = vtail_horses
                         )

wing_mac = mean_aerodynamic_center(wing)
x_w      = wing_mac[1]

## Case
ac_name = "My Aircraft"
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
ρ       = 1.225
ref     = [ x_w, 0., 0.]
V, α, β = 1.0, 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

@time data =
    solve_case(aircraft, fs;
               rho_ref     = ρ,         # Reference density
               r_ref       = ref,       # Reference point for moments
               area_ref    = S,         # Reference area
               span_ref    = b,         # Reference span
               chord_ref   = c,         # Reference chord
               name        = ac_name,   # Aircraft name
               print       = true,      # Prints the results for the entire aircraft
               print_components = true, # Prints the results for each component
              );

## Data access
Γs_comp    = data.wing.circulations
CFs        = data.wing.CFs
CMs        = data.wing.CMs

# Component data processing
nf_coeffs  = nearfield(data.wing) # Nearfield coefficients
ff_coeffs  = farfield(data.wing)  # Farfield coefficients

print_coefficients(nf_coeffs, ff_coeffs, :wing)

## System data processing
Γs   = circulations(data)                # Circulations of entire system
nf_c = sum(nearfield_coefficients(data)) # Total nearfield force and moment coefficients
ff_c = sum(farfield_coefficients(data))  # Total farfield force coefficients

print_coefficients(nf_c, ff_c, ac_name)

## Stability case
@time dv_data =
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

## Streamlines

# Chordwise distribution
# using StaticArrays

# num_points = 100
# max_z = 2
# y = span(wing) / 10
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 10
init        = chop_leading_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz]) ;
                init .+ Ref([dx, dy,-dz]) ];

distance = 8
num_stream_points = 200
streams = plot_streams(fs, seed, aircraft, Γs, distance, num_stream_points);

##
using Plots
gr(size = (1280, 720), dpi = 300)

##
horseshoe_coords = plot_panels(aircraft_panels)
horseshoe_points = Tuple.(horseshoe_point.(aircraft_panels));

##
z_limit = b
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     camera = (30, 60),
    #  xlim = (-z_limit/2, z_limit/2),
     aspect_ratio = 1,
     zlim = (-z_limit/2, z_limit/2),
     size = (1280, 720)
    )
plot!.(horseshoe_coords, color = :gray, label = :none)
scatter!(horseshoe_points, marker = 1, color = :black, label = :none)
plot!.(streams, color = :green, label = :none)
plot!()

## Exaggerated CF distribution for plot, only works with GR and not Plotly

# Forces
wind_CFs = body_to_wind_axes.(CFs, fs.alpha, fs.beta)
CDis     = @. getindex(wind_CFs, 1)
CYs      = @. getindex(wind_CFs, 2)
CLs      = @. getindex(wind_CFs, 3)

hs_pts = Tuple.(bound_leg_center.(aircraft))

quiver!(hs_pts, quiver=(CDis[:], CYs[:], CLs[:]) .* 500)
plot!()