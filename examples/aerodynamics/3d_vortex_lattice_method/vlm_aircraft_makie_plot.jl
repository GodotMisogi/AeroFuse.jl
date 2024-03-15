## NOTE: This is called in vlm_aircraft.jl. It won't work by itself.
using CairoMakie
CairoMakie.activate!()

# using WGLMakie
# WGLMakie.activate!()

set_theme!(
# theme_black()
# theme_light()
)

## Streamlines
# Spanwise distribution
span_points = 20
dx, dy, dz = 0, 0, 1e-3
init = chop_leading_edge(wing, span_points)
seed = init .+ Ref([dx, dy, dz])
#   init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = streamlines(sys, seed, distance, num_stream_points);

wing_cam_connec = triangle_connectivities(LinearIndices(camber_coordinates(wing_mesh)))
htail_cam_connec = triangle_connectivities(LinearIndices(camber_coordinates(htail_mesh)))
vtail_cam_connec = triangle_connectivities(LinearIndices(camber_coordinates(vtail_mesh)));

## Surface velocities
vels = surface_velocities(sys);
sps = norm.(vels)

wing_sp_points = extrapolate_point_mesh(sps.wing)
htail_sp_points = extrapolate_point_mesh(sps.htail)
vtail_sp_points = extrapolate_point_mesh(sps.vtail)

## Surface pressure coefficients
cps = norm.(CFs) * S

wing_cp_points = extrapolate_point_mesh(cps.wing)
htail_cp_points = extrapolate_point_mesh(cps.htail)
vtail_cp_points = extrapolate_point_mesh(cps.vtail)

## Figure plot
fig1 = Figure(resolution = (1280, 720))

scene = LScene(fig1[1:4, 1])
ax = fig1[1:4, 2] = GridLayout()

ax_cd = Axis(ax[1, 1:2], ylabel = L"C_{D_i}", title = LS("Spanwise Loading"))
ax_cy = Axis(ax[2, 1:2], ylabel = L"C_Y")
ax_cl = Axis(ax[3, 1:2], xlabel = L"y", ylabel = L"C_L")

# Spanload plot
function plot_spanload!(ax, ll_loads, name = "Wing")
    @views lines!(ax[1, 1:2], ll_loads[:, 1], ll_loads[:, 2], label = name)
    @views lines!(ax[2, 1:2], ll_loads[:, 1], ll_loads[:, 3], label = name)
    @views lines!(ax[3, 1:2], ll_loads[:, 1], ll_loads[:, 4], label = name)

    nothing
end

# Surface pressure meshes
m1 = poly!(
    scene,
    vec(camber_coordinates(wing_mesh)),
    wing_cam_connec,
    color = vec(wing_cp_points),
    colormap = :rainbow1,
)
m2 = poly!(
    scene,
    vec(camber_coordinates(htail_mesh)),
    htail_cam_connec,
    color = vec(htail_cp_points),
    colormap = :rainbow1,
)
m3 = poly!(
    scene,
    vec(camber_coordinates(vtail_mesh)),
    vtail_cam_connec,
    color = vec(vtail_cp_points),
    colormap = :rainbow1,
)

Colorbar(fig1[4, 1], m1.plots[1], label = L"Pressure Coefficient, $C_p$", vertical = false)

fig1[0, :] = Label(fig1, LS("Vortex Lattice Analysis"), fontsize = 20)

# Airfoil meshes
# wing_surf = surface_coordinates(wing_mesh, wing_mesh.n_span, 60)
# surf_connec = triangle_connectivities(LinearIndices(wing_surf))
# wing_surf_mesh = mesh(vec(wing_surf), surf_connec)
# w1 = wireframe!(scene, wing_surf_mesh.plot[1][], color = :grey, alpha = 0.1)

# Borders
lines!(scene, plot_planform(wing))
lines!(scene, plot_planform(htail))
lines!(scene, plot_planform(vtail))

# Streamlines
[lines!(scene, Point3f.(stream[:]), color = :lightblue) for stream in eachcol(streams)]

plot_spanload!(ax, wing_ll, LS("Wing"))
plot_spanload!(ax, htail_ll, LS("Horizontal Tail"))
plot_spanload!(ax, vtail_ll, LS("Vertical Tail"))

# Legend
Legend(ax[4, 1:2], ax_cl)

fig1.scene

## Save figure
# save("plots/VortexLattice.pdf", fig1, px_per_unit = 1.5)

## Animation
#=======================================================#

# pts = [ Point3f[stream] for stream in streams[1,:] ]

# [ lines!(scene, pts[i], color = :green) for i in eachindex(pts) ]

# Recording (Doesn't work right now)
# fps     = 30
# nframes = length(streams[:,1])

# record(fig1, "plots/vlm_animation.mp4", 1:nframes) do i 
#     for j in eachindex(streams[1,:])
#         pts[j][] = push!(pts[j][], Point3f(streams[i,j]))
#     end
#     sleep(1/fps) # refreshes the display!
#     notify(pts[i])
# end