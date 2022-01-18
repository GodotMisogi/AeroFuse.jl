## Wing analysis case
using AeroMDAO

## Wing section setup
wing_right = HalfWing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
                      chords    = [1.0, 0.6, 0.2],
                      twists    = [0.0, 0.0, 0.0],
                      spans     = [5.0, 0.5],
                      dihedrals = [5., 5.],
                      LE_sweeps = [5., 5.]);
wing = Wing(wing_right, wing_right)
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Meshing and assembly
wing_mesh    = WingMesh(wing, [12, 3], 6, 
                        span_spacing = Cosine()
                       );
aircraft     = ComponentVector(wing = make_horseshoes(wing_mesh))

# Freestream conditions
fs      = Freestream(alpha = 1.0, 
                     beta  = 0.0, 
                     omega = [0.,0.,0.])

# Reference values
refs    = References(
                     speed     = 10.0,
                     density   = 1.225,
                     viscosity = 1.5e-5,
                     area      = projected_area(wing),
                     span      = span(wing), 
                     chord     = mean_aerodynamic_chord(wing), 
                     location  = mean_aerodynamic_center(wing)
                    )

## Solve system
system  = solve_case(aircraft, fs, refs; 
                     print = true
                    )

## Compute dynamics
ax       = Wind()
CFs, CMs = surface_coefficients(system; axes = ax)
# Fs       = surface_forces(system)
# Ms       = surface_moments(system)
# Fs, Ms   = surface_dynamics(system; axes = ax) 

## Viscous drag prediction using empirical models

## Equivalent flat-plate skin-friction estimation
x_tr        = [0.35, 0.35]              # Transition locations over sections
CDv_plate   = profile_drag_coefficient(wing, x_tr, refs.speed, refs.density, 330., refs.area, refs.viscosity)

## Local-dissipation drag estimation (WRONG???)
cam_panels  = camber_panels(wing_mesh)
edge_speeds = surface_velocities(system).wing
M           = mach_number(system.reference)  # Mach number
CDv_diss    = local_dissipation_drag(wing, panel_area.(cam_panels), refs.density, edge_speeds, [0.2, 0.35, 0.35, 0.2], refs.speed, refs.density, M, refs.viscosity) / system.reference.area

## Viscous drag coefficient
CDv = CDv_plate

## Total force coefficients with viscous drag prediction
CDi_nf, CY_nf, CL_nf, Cl, Cm, Cn = nf = nearfield(system) 
CDi_ff, CY_ff, CL_ff = ff = farfield(system)

nf_v = [ CDi_nf + CDv; CDv; nf ]
ff_v = [ CDi_ff + CDv; CDv; ff ]
print_coefficients(nf_v, ff_v, :wing)


## Evaluate case with stability derivatives
@time dv_data = solve_case_derivatives(aircraft, fs, refs;
                                       print = true
                                      );

## Plotting
# using StaticArrays
using Plots
gr()

## Streamlines

# Chordwise distribution
# num_points = 50
# max_z = 0.1
# y = span(wing) / 2 - 0.05
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

# Spanwise distribution
span_points = 20
init        = chop_trailing_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy,  dz])  ;
                init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = plot_streamlines(system, seed, distance, num_stream_points);

## Display
horseshoe_panels = chord_panels(wing_mesh)
horseshoe_coords = plot_panels(horseshoe_panels)
wing_coords      = plot_wing(wing);
horseshoe_points = Tuple.(horseshoe_point.(system.vortices))
ys               = getindex.(horseshoe_points, 2)

## Coordinates
z_limit = 5
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (30, 60),
     zlim = (-0.1, z_limit),
     size = (800, 600))
# plot!(wing_coords[:,1], wing_coords[:,2], wing_coords[:,3])
plot!.(horseshoe_coords, color = :black, label = :none)
scatter!(vec(horseshoe_points), marker = 1, color = :black, label = :none)
[ plot!(stream, color = :green, label = :none) for stream in eachcol(Tuple.(streams)) ]
plot!()

## Spanwise forces
LL_loads    = span_loads(horseshoe_panels, CFs.wing, projected_area(wing))
CL_loadings = vec(sum(system.circulations.wing, dims = 1)) / (0.5 * refs.speed * refs.chord)

## Lifting line loads
plot_CD = plot(LL_loads[:,1], LL_loads[:,2], label = :none, ylabel = "CDi")
plot_CY = plot(LL_loads[:,1], LL_loads[:,3], label = :none, ylabel = "CY")
plot_CL = begin
            plot(LL_loads[:,1], LL_loads[:,4], label = :none, xlabel = "y", ylabel = "CL")
            plot!(LL_loads[:,1], CL_loadings, label = "Normalized", xlabel = "y")
          end
plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))

## Lift distribution

# Exaggerated CF distribution for plot
hs_pts = vec(Tuple.(bound_leg_center.(system.vortices)))

plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (60, 60),
     zlim = (-0.1, z_limit),
     title = "Forces (Exaggerated)"
    )
plot!.(horseshoe_coords, color = :gray, label = :none)
# scatter!(cz_pts, zcolor = vec(CLs), marker = 2, label = "CL (Exaggerated)")
quiver!(hs_pts, quiver=(LL_loads[:,2], LL_loads[:,3], LL_loads[:,4]) .* 10)
plot!(size = (800, 600))