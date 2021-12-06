## Wing analysis case
using Revise
using AeroMDAO
using ComponentArrays
using LinearAlgebra
using SplitApplyCombine

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
                          wing  = make_vortex_rings(wing_mesh),
                          htail = make_vortex_rings(htail_mesh),
                          vtail = make_vortex_rings(vtail_mesh)
                         );

## Case
ac_name = :aircraft
ρ       = 1.225
ref     = [ x_w, 0., 0.]
V, α, β = 15.0, 0.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)
refs    = References(S, b, c, ρ, ref)

##
@time begin 
    data = solve_case(aircraft, fs, refs;
                      print            = true, # Prints the results for only the aircraft
                      print_components = true, # Prints the results for all components
                    #   finite_core      = true
                     );

    # Compute dynamics
    ax       = Stability()
    Fs       = surface_forces(data)
    Fs, Ms   = surface_dynamics(data; axes = ax) 
    CFs, CMs = surface_coefficients(data; axes = ax)

    nfs = nearfield_coefficients(data)
    ffs = farfield_coefficients(data)

    nf  = nearfield(data) 
    ff  = farfield(data)
end;

## Spanwise forces
using StaticArrays
function lifting_line_loads(panels, CFs, S)
    CFs  = combinedimsview(CFs)
    CDis = @views CFs[1,:,:]
    CYs  = @views CFs[2,:,:]
    CLs  = @views CFs[3,:,:]

    area_scale  = S ./ sum(panel_area, panels, dims = 1)[:]
    span_CDis   = sum(CDis, dims = 1)[:] .* area_scale
    span_CYs    = sum(CYs,  dims = 1)[:] .* area_scale
    span_CLs    = sum(CLs,  dims = 1)[:] .* area_scale

    # ys = sum(x -> getindex(midpoint(x), 2), panels, dims = 1) 

    SVector.(span_CDis, span_CYs, span_CLs)
end

hs_pts   = horseshoe_point.(data.horseshoes)
wing_ys  = @views getindex.(hs_pts.wing[1,:],  2)
htail_ys = @views getindex.(hs_pts.htail[1,:], 2)
vtail_ys = @views getindex.(hs_pts.vtail[1,:], 2)

##
wing_ll  = lifting_line_loads(chord_panels(wing_mesh), CFs.wing, S)
htail_ll = lifting_line_loads(chord_panels(htail_mesh), CFs.htail, S)
vtail_ll = lifting_line_loads(chord_panels(vtail_mesh), CFs.vtail, S);

## Plotting
using GLMakie
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
streams = plot_streams(fs, seed, data.horseshoes, data.circulations, distance, num_stream_points);

## Mesh connectivities
triangle_connectivities(inds) = @views [ inds[1:end-1,1:end-1][:] inds[1:end-1,2:end][:]   inds[2:end,2:end][:]   ;
                                           inds[2:end,2:end][:]   inds[2:end,1:end-1][:] inds[1:end-1,1:end-1][:] ]

wing_cam_connec  = triangle_connectivities(LinearIndices(wing_mesh.cam_mesh))
htail_cam_connec = triangle_connectivities(LinearIndices(htail_mesh.cam_mesh))
vtail_cam_connec = triangle_connectivities(LinearIndices(vtail_mesh.cam_mesh));

## Extrapolating surface values to neighbouring points
function extrapolate_point_mesh(mesh)
    m, n   = size(mesh)
    points = zeros(eltype(mesh), m + 1, n + 1)

    # The quantities are measured at the bound leg (0.25×)
    @views points[1:end-1,1:end-1] += 0.75 * mesh / 2
    @views points[1:end-1,2:end]   += 0.75 * mesh / 2
    @views points[2:end,1:end-1]   += 0.25 * mesh / 2
    @views points[2:end,2:end]     += 0.25 * mesh / 2

    points
end

## Surface velocities
vels = surface_velocities(data);
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
fig  = Figure(resolution = (1280, 720))

scene = LScene(fig[1:4,1])
ax    = fig[1:4,2] = GridLayout()

ax_cd = GLMakie.Axis(ax[1,1:2], ylabel = L"C_{D_i}", title = LS("Spanwise Loading"))
ax_cy = GLMakie.Axis(ax[2,1:2], ylabel = L"C_Y",)
ax_cl = GLMakie.Axis(ax[3,1:2], xlabel = L"y", ylabel = L"C_L")

# Spanload plot
function plot_spanload!(ax, ys, CFs, name = "Wing")
    CFs  = combinedimsview(CFs)
    @views lines!(ax[1,1:2], ys, CFs[1,:], label = name,)
    @views lines!(ax[2,1:2], ys, CFs[2,:], label = name,)
    @views lines!(ax[3,1:2], ys, CFs[3,:], label = name,)

    nothing
end

plot_spanload!(ax, wing_ys, wing_ll, LS("Wing"))
plot_spanload!(ax, htail_ys, htail_ll, LS("Horizontal Tail"))
plot_spanload!(ax, vtail_ys, vtail_ll, LS("Vertical Tail"))

# Legend
Legend(ax[4,1:2], ax_cl)
fig[0, :] = Label(fig, LS("Vortex Lattice Analysis"), textsize = 20)

# Surface pressure meshes
m1 = poly!(scene, wing_mesh.cam_mesh[:],  wing_cam_connec,  color =  wing_cp_points[:])
m2 = poly!(scene, htail_mesh.cam_mesh[:], htail_cam_connec, color = htail_cp_points[:])
m3 = poly!(scene, vtail_mesh.cam_mesh[:], vtail_cam_connec, color = vtail_cp_points[:])

# Airfoil meshes
# wing_surf = surface_coordinates(wing_mesh, wing_mesh.n_span, 60)
# surf_connec = triangle_connectivities(LinearIndices(wing_surf))
# wing_surf_mesh = mesh(wing_surf[:], surf_connec)
# w1 = wireframe!(scene, wing_surf_mesh.plot[1][], color = :grey, alpha = 0.1)

# Borders
lines!(scene, plot_wing(wing))
lines!(scene, plot_wing(htail))
lines!(scene, plot_wing(vtail))

l1 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(make_panels(wing_mesh.cam_mesh))  ]
l2 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(make_panels(htail_mesh.cam_mesh)) ]
l3 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(make_panels(vtail_mesh.cam_mesh)) ]

# Streamlines
[ lines!(scene, stream[:], color = :green) for stream in eachcol(streams) ]

fig.scene

## Save figure
# save("plots/VortexLattice.png", fig, px_per_unit = 1.5)

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
# hs_pts = Tuple.(bound_leg_center.(horses))[:]
# arrows!(scene, getindex.(hs_pts, 1), getindex.(hs_pts, 2), getindex.(hs_pts, 3), 
#                 CDis[:], CYs[:], CLs[:], 
#                 arrowsize = Vec3f.(0.3, 0.3, 0.4),
#                 lengthscale = 10,
#                 label = "Forces (Exaggerated)")