## Wing analysis case
using AeroMDAO
import LinearAlgebra: norm

## Surfaces

# Wing
wing = Wing(foils     = Foil.(fill(naca4(2,4,1,2), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            LE_sweeps = [5.]);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4(0,0,1,2), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4(0,0,0,9), 2)),
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

## WingMesh type
wing_mesh  = WingMesh(wing, [12], 6, 
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
                     area     = projected_area(wing),
                     span     = span(wing),
                     chord    = mean_aerodynamic_chord(wing),
                     density  = 1.225,
                     location = [ x_w, 0., 0.])

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

## Spanwise forces/lifting line loads
wing_ll  = span_loads(chord_panels(wing_mesh), CFs.wing, S)
htail_ll = span_loads(chord_panels(htail_mesh), CFs.htail, S)
vtail_ll = span_loads(chord_panels(vtail_mesh), CFs.vtail, S);

## Plotting
using CairoMakie
using LaTeXStrings

set_theme!(
            # theme_black()
            # theme_light()
          )

const LS = LaTeXString

## Streamlines
# Spanwise distribution
span_points = 30
init        = chop_leading_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-1
seed        = init # [ init .+ Ref([dx, dy,  dz])
                   #   init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = streamlines(system, seed, distance, num_stream_points);

wing_cam_connec  = triangle_connectivities(LinearIndices(wing_mesh.camber_mesh))
htail_cam_connec = triangle_connectivities(LinearIndices(htail_mesh.camber_mesh))
vtail_cam_connec = triangle_connectivities(LinearIndices(vtail_mesh.camber_mesh));

## Surface velocities
vels = surface_velocities(system);
sps  = norm.(vels)

wing_sp_points  = extrapolate_point_mesh(sps.wing)
htail_sp_points = extrapolate_point_mesh(sps.htail)
vtail_sp_points = extrapolate_point_mesh(sps.vtail)

## Surface pressure coefficients
cps  = norm.(CFs) * S

wing_cp_points  = extrapolate_point_mesh(cps.wing)
htail_cp_points = extrapolate_point_mesh(cps.htail)
vtail_cp_points = extrapolate_point_mesh(cps.vtail)

## Figure plot
fig1  = Figure(resolution = (1280, 720))

scene = LScene(fig1[1:4,1])
ax    = fig1[1:4,2] = GridLayout()

ax_cd = Axis(ax[1,1:2], ylabel = L"C_{D_i}", title = LS("Spanwise Loading"))
ax_cy = Axis(ax[2,1:2], ylabel = L"C_Y",)
ax_cl = Axis(ax[3,1:2], xlabel = L"y", ylabel = L"C_L")

# Spanload plot
function plot_spanload!(ax, ll_loads, name = "Wing")
    @views lines!(ax[1,1:2], ll_loads[:,1], ll_loads[:,2], label = name,)
    @views lines!(ax[2,1:2], ll_loads[:,1], ll_loads[:,3], label = name,)
    @views lines!(ax[3,1:2], ll_loads[:,1], ll_loads[:,4], label = name,)

    nothing
end

plot_spanload!(ax, wing_ll, LS("Wing"))
plot_spanload!(ax, htail_ll, LS("Horizontal Tail"))
plot_spanload!(ax, vtail_ll, LS("Vertical Tail"))

# Legend
Legend(ax[4,1:2], ax_cl)
fig1[0, :] = Label(fig1, LS("Vortex Lattice Analysis"), textsize = 20)

# Surface pressure meshes
m1 = poly!(scene, vec(wing_mesh.camber_mesh),  wing_cam_connec,  color = vec(wing_cp_points))
m2 = poly!(scene, vec(htail_mesh.camber_mesh), htail_cam_connec, color = vec(htail_cp_points))
m3 = poly!(scene, vec(vtail_mesh.camber_mesh), vtail_cam_connec, color = vec(vtail_cp_points))

# Airfoil meshes
# wing_surf = surface_coordinates(wing_mesh, wing_mesh.n_span, 60)
# surf_connec = triangle_connectivities(LinearIndices(wing_surf))
# wing_surf_mesh = mesh(vec(wing_surf), surf_connec)
# w1 = wireframe!(scene, wing_surf_mesh.plot[1][], color = :grey, alpha = 0.1)

# Borders
lines!(scene, plot_wing(wing))
lines!(scene, plot_wing(htail))
lines!(scene, plot_wing(vtail))

# l1 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(wing_mesh))  ]
# l2 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(htail_mesh)) ]
# l3 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(vtail_mesh)) ]

# Streamlines
[ lines!(scene, stream[:], color = :green) for stream in eachcol(streams) ]

fig1.scene

## Save figure
save("plots/VortexLattice.svg", fig1, px_per_unit = 1.5)

## Animation settings
pts = [ Node(Point3f0[stream]) for stream in streams[1,:] ]

[ lines!(scene, pts[i], color = :green, axis = (; type = Axis3)) for i in eachindex(pts) ]

# Recording
fps     = 30
nframes = length(streams[:,1])

record(fig, "plots/vlm_animation.mp4", 1:nframes) do i 
    for j in eachindex(streams[1,:])
        pts[j][] = push!(pts[j][], Point3f0(streams[i,j]))
    end
    sleep(1/fps) # refreshes the display!
    notify(pts[i])
end

## Arrows
# hs_pts = vec(Tuple.(bound_leg_center.(horses)))
# arrows!(scene, getindex.(hs_pts, 1), getindex.(hs_pts, 2), getindex.(hs_pts, 3), 
#                 vec(CDis),vec(CYs),vec( Ls), 
#                 arrowsize = Vec3f.(0.3, 0.3, 0.4),
#                 lengthscale = 10,
#                 label = "Forces (Exaggerated)")