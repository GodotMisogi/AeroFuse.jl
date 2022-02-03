## Aircraft analysis case
using AeroMDAO

## Lifting surfaces

# Wing
wing = Wing(foils     = fill(naca4((0,0,1,2)), 3),
            chords    = [1.0, 0.6, 0.2],
            twists    = [2.0, 0.0, -2.0],
            spans     = [4.0, 0.2],
            dihedrals = [5., 30.],
            sweeps      = [5., 30.]);

# Horizontal tail
htail = Wing(foils     = fill(naca4(0,0,1,2), 2),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweeps      = [6.39],
             position  = [4., 0, 0.2],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = fill(naca4(0,0,0,9), 2),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweeps      = [7.97],
                 position  = [4., 0, 0.2],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## WingMesh type
wing_mesh  = WingMesh(wing, [12, 4], 6, 
                      span_spacing = Cosine()
                     )
htail_mesh = WingMesh(htail, [12], 6, 
                      span_spacing = Cosine()
                     )
vtail_mesh = WingMesh(vtail, [12], 6, 
                      span_spacing = Cosine()
                     )

aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );

## Case
fs      = Freestream(alpha = 0.0, 
                     beta  = 0.0, 
                     omega = [0., 0., 0.]);

refs    = References(speed    = 1.0, 
                     density  = 1.225,
                     area     = projected_area(wing),
                     span     = span(wing),
                     chord    = mean_aerodynamic_chord(wing),
                     location = mean_aerodynamic_center(wing))

##
@time begin 
    system = solve_case(aircraft, fs, refs;
                        print            = true, # Prints the results for only the aircraft
                        print_components = true, # Prints the results for all components
                      #   finite_core      = true
                       );

    # Compute dynamics
    ax       = Geometry() # Geometry, Stability(), Body()
    CFs, CMs = surface_coefficients(system; axes = ax)
    Fs, Ms   = surface_dynamics(system; axes = ax)
    # Fs       = surface_forces(system; axes = ax)
    # vels     = surface_velocities(system)

    nfs = nearfield_coefficients(system)
    ffs = farfield_coefficients(system)

    nf  = nearfield(system) 
    ff  = farfield(system)
end;

## Streamlines

# Spanwise distribution
span_points = 10
init        = chop_leading_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz]) ;
                init .+ Ref([dx, dy,-dz]) ];

distance = 5
num_stream_points = 100
streams = streamlines(system, seed, distance, num_stream_points);

## Fuselage definition
lens = [0.0, 0.05,0.02,  0.3, 0.6, 0.8, 1.0]
rads = [0.0, 0.3, 0.25, 0.4, 0.2, 0.1, 0.00]
fuse = Fuselage(5.5, lens, rads, [-1, 0.,0.2])

lens_rads = coordinates(fuse, 40)

## Circles for plotting
n_pts          = 20
circle3D(r, n) = let arcs = 0:2π/n:2π; [ zeros(length(arcs)) r * cos.(arcs) r * sin.(arcs) ] end

xs = lens_rads[:,1]
circs = [ reduce(hcat, eachrow(circ) .+ Ref([x; 0; 0] + fuse.position))' for (x, circ) in zip(xs, circle3D.(lens_rads[:,2], n_pts)) ]

##
using Plots
gr(size = (1280, 720), dpi = 300)
# plotlyjs(dpi = 150)

##
aircraft_panels   = ComponentArray(wing  = chord_panels(wing_mesh),
                                   htail = chord_panels(htail_mesh),
                                   vtail = chord_panels(vtail_mesh))
panel_coordinates = plot_panels(aircraft_panels)
horseshoe_points  = Tuple.(horseshoe_point.(aircraft_panels))[:];

wing_coords  = plot_wing(wing)
htail_coords = plot_wing(htail)
vtail_coords = plot_wing(vtail)

z_limit = span(wing)
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     camera = (30, 60),
    #  xlim = (-z_limit/2, z_limit/2),
    #  aspect_ratio = 1,
     zlim = (-z_limit/2, z_limit/2),
     size = (1280, 720)
    )
plot!.(panel_coordinates, color = :gray, label = :none)
plot!(wing_coords, label = "Wing")
plot!(htail_coords, label = "Horizontal Tail")
plot!(vtail_coords, label = "Vertical Tail")
# scatter!(horseshoe_points, marker = 1, color = :black, label = :none)
[ plot!(Tuple.(stream), color = :green, label = :none) for stream in eachcol(streams) ]
plot!()

[ plot!(circ[:,1], circ[:,2], circ[:,3], color = :gray, label = :none) for circ in circs ] 
plot!()

## Exaggerated CF distribution for plot, only works with GR and not Plotly

# Forces
wind_CFs = geometry_to_wind_axes.(CFs, fs.alpha, fs.beta)
CDis     = @. getindex(wind_CFs, 1)
CYs      = @. getindex(wind_CFs, 2)
CLs      = @. getindex(wind_CFs, 3)

hs_pts = Tuple.(bound_leg_center.(horses))[:]

quiver!(hs_pts, quiver=(CDis[:], CYs[:], CLs[:]) .* 500)
plot!()