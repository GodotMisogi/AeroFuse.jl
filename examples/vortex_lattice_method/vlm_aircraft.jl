## Wing analysis case
using AeroMDAO

## Surfaces

# Wing
wing = Wing(
    foils     = fill(naca4(2,4,1,2), 2),
    chords    = [1.0, 0.6],
    twists    = [2.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.],
    sweeps    = [5.],
    symmetry  = true,
    # flip      = true
);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25],
    dihedrals = [0.],
    sweeps    = [6.39],
    position  = [4., 0, -0.1],
    angle     = -2.,
    axis      = [0., 1., 0.],
    symmetry  = true
)

# Vertical tail
vtail = Wing(
    foils     = fill(naca4(0,0,0,9), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.0],
    dihedrals = [0.],
    sweeps    = [7.97],
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.]
)

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## WingMesh type
wing_mesh  = WingMesh(wing, [24], 6,
    # span_spacing = Cosine()
)
htail_mesh = WingMesh(htail, [24], 6, 
    # span_spacing = Cosine()
)
vtail_mesh = WingMesh(vtail, [24], 6, 
    span_spacing = Cosine()
)

aircraft =  ComponentVector(
    wing  = make_horseshoes(wing_mesh),
    htail = make_horseshoes(htail_mesh),
    vtail = make_horseshoes(vtail_mesh)
);

## Case
fs  = Freestream(
    alpha = 0.0, 
    beta  = 0.0, 
    omega = [0., 0., 0.]
);

ref = References(
    speed     = 1.0, 
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing)
)

##
@time begin 
    system = solve_case(aircraft, fs, ref;
                        print            = true, # Prints the results for only the aircraft
                        print_components = true, # Prints the results for all components
                       );

    # Compute dynamics
    ax       = Geometry() # Geometry, Stability(), Body()
    CFs, CMs = surface_coefficients(system; axes = ax)
    Fs, Ms   = surface_dynamics(system; axes = ax)
    # Fs       = surface_forces(system; axes = ax)
    # vels     = surface_velocities(system)

    nfs = nearfield_coefficients(system)
    ffs = farfield_coefficients(system)
end;

## Viscous drag prediction

# Equivalent flat-plate skin friction estimation
CDv_wing  = profile_drag_coefficient(wing,  [0.8, 0.8], system.reference)
CDv_htail = profile_drag_coefficient(htail, [0.6, 0.6], system.reference)
CDv_vtail = profile_drag_coefficient(vtail, [0.6, 0.6], system.reference)

CDv_plate = CDv_wing + CDv_htail + CDv_vtail

## Local dissipation form factor friction estimation (WRONG???)
import LinearAlgebra: norm

edge_speeds = norm.(surface_velocities(system)); # Inviscid speeds on the surfaces

# Drag coefficients
CDvd_wing  = profile_drag_coefficient(wing_mesh,  [0.8, 0.8], edge_speeds.wing,  ref)
CDvd_htail = profile_drag_coefficient(htail_mesh, [0.6, 0.6], edge_speeds.htail, ref)
CDvd_vtail = profile_drag_coefficient(vtail_mesh, [0.6, 0.6], edge_speeds.vtail, ref)

CDv_diss = CDvd_wing + CDvd_htail + CDvd_vtail

## Viscous drag coefficient
CDv = CDv_plate

## Total force coefficients with empirical viscous drag prediction
CDi_nf, CY_nf, CL_nf, Cl, Cm, Cn = nf = nearfield(system) 
CDi_ff, CY_ff, CL_ff = ff = farfield(system)

nf_v = (CD = CDi_nf + CDv, CDv = CDv, nf...)
ff_v = (CD = CDi_ff + CDv, CDv = CDv, ff...)

## Spanwise forces/lifting line loads
wing_ll  = spanwise_loading(wing_mesh, CFs.wing,  S)
htail_ll = spanwise_loading(htail_mesh, CFs.htail, S)
vtail_ll = spanwise_loading(vtail_mesh, CFs.vtail, S);

## Plotting
using CairoMakie
CairoMakie.activate!()

set_theme!(
            # theme_black()
            # theme_light()
          )

using LaTeXStrings
const LS = LaTeXString

## Streamlines
# Spanwise distribution
span_points = 20
init        = chop_leading_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = init .+ Ref([dx, dy,  dz])
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

# Surface pressure meshes
m1 = poly!(scene, vec(wing_mesh.camber_mesh),  wing_cam_connec,  color = vec(wing_cp_points))
m2 = poly!(scene, vec(htail_mesh.camber_mesh), htail_cam_connec, color = vec(htail_cp_points))
m3 = poly!(scene, vec(vtail_mesh.camber_mesh), vtail_cam_connec, color = vec(vtail_cp_points))

Colorbar(fig1[4,1], m1.plots[1], label = L"Pressure Coefficient, $C_p$", vertical = false)

fig1[0, :] = Label(fig1, LS("Vortex Lattice Analysis"), textsize = 20)


# Airfoil meshes
# wing_surf = surface_coordinates(wing_mesh, wing_mesh.n_span, 60)
# surf_connec = triangle_connectivities(LinearIndices(wing_surf))
# wing_surf_mesh = mesh(vec(wing_surf), surf_connec)
# w1 = wireframe!(scene, wing_surf_mesh.plot[1][], color = :grey, alpha = 0.1)

# Borders
lines!(scene, plot_planform(wing))
lines!(scene, plot_planform(htail))
lines!(scene, plot_planform(vtail))

# l1 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(wing_mesh))  ]
# l2 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(htail_mesh)) ]
# l3 = [ lines!(scene, pts, color = :grey) for pts in plot_panels(camber_panels(vtail_mesh)) ]

# Streamlines
[ lines!(scene, Point3f.(stream[:]), color = :green) for stream in eachcol(streams) ]

fig1.scene

## Save figure
# save("plots/VortexLattice.pdf", fig1, px_per_unit = 1.5)

## Animation settings
pts = [ Node(Point3f0[stream]) for stream in streams[1,:] ]

[ lines!(scene, pts[i], color = :green, axis = (; type = Axis3)) for i in eachindex(pts) ]

# Recording
fps     = 30
nframes = length(streams[:,1])

# record(fig, "plots/vlm_animation.mp4", 1:nframes) do i 
#     for j in eachindex(streams[1,:])
#         pts[j][] = push!(pts[j][], Point3f0(streams[i,j]))
#     end
#     sleep(1/fps) # refreshes the display!
#     notify(pts[i])
# end

## Arrows
# hs_pts = vec(Tuple.(bound_leg_center.(horses)))
# arrows!(scene, getindex.(hs_pts, 1), getindex.(hs_pts, 2), getindex.(hs_pts, 3), 
#                 vec(CDis),vec(CYs),vec( Ls), 
#                 arrowsize = Vec3f.(0.3, 0.3, 0.4),
#                 lengthscale = 10,
#                 label = "Forces (Exaggerated)")

## VARIABLE ANALYSES
#=========================================================#

using Setfield
using Base.Iterators: product
using Plots

pyplot()

## Speed sweep
Vs = 1.0:10:300
res_Vs = combinedimsview(
    map(Vs) do V
        ref1 = @set ref.speed = V
        sys = solve_case(aircraft, fs, ref1)
        [ V; farfield(sys)...; nearfield(sys)... ]
    end, (1)
)

## 
Plots.plot(
    res_Vs[:,1], res_Vs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "V",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)

##
f2 = Figure()
ax1 = Axis(f2[1,1], xlabel = L"V", ylabel = L"C_{D_{ff}}")
lines!(res_Vs[:,1], res_Vs[:,2]) # CD_ff
ax2 = Axis(f2[1,2], xlabel = L"V", ylabel = L"C_{Y_{ff}}")
lines!(res_Vs[:,1], res_Vs[:,3])
ax3 = Axis(f2[1,3], xlabel = L"V", ylabel = L"CL_{ff}")
lines!(res_Vs[:,1], res_Vs[:,4])
ax4 = Axis(f2[2,1], xlabel = L"V", ylabel = L"C_D")
lines!(res_Vs[:,1], res_Vs[:,5])
ax5 = Axis(f2[2,2], xlabel = L"V", ylabel = L"C_Y")
lines!(res_Vs[:,1], res_Vs[:,6])
ax6 = Axis(f2[2,3], xlabel = L"V", ylabel = L"C_L")
lines!(res_Vs[:,1], res_Vs[:,7])
ax7 = Axis(f2[3,1], xlabel = L"V", ylabel = L"C_\ell")
lines!(res_Vs[:,1], res_Vs[:,8])
ax8 = Axis(f2[3,2], xlabel = L"V", ylabel = L"C_m")
lines!(res_Vs[:,1], res_Vs[:,9])
ax9 = Axis(f2[3,3], xlabel = L"V", ylabel = L"C_n")
lines!(res_Vs[:,1], res_Vs[:,10])

f2

## Alpha sweep
αs = -5:0.5:5
res_αs = combinedimsview(
    map(αs) do α
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref)
        [ α; farfield(sys)...; nearfield(sys)... ]
    end, (1)
)

Plots.plot(
    res_αs[:,1], res_αs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "α",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)

## (Speed, alpha) sweep
res = combinedimsview(
    map(product(Vs, αs)) do (V, α)
        ref1 = @set ref.speed = V
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref1)
        [ V; α; farfield(sys)...; nearfield(sys)... ]
    end
)

##
res_p = permutedims(res, (3,1,2))

# CDi
plt_CDi_ff = Plots.plot(camera = (75, 30), ylabel = L"\alpha", xlabel = "V", zlabel = L"C_{D_i}")
[ Plots.plot!(res_p[:,1,n], res_p[:,2,n], res_p[:,3,n], label = "", c = :black) for n in axes(res_p,3) ]

# CL
plt_CL_ff = Plots.plot(camera = (75,15), 
ylabel = L"\alpha", xlabel = "V", zlabel = L"C_{L}", 
label = "", c = :black)
[ Plots.plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,8,n], label = "", c = :black) for n in axes(res_p,3) ]

plt_Cm_ff = Plots.plot(camera = (100,30), ylabel = L"\alpha", xlabel = "V", zlabel = L"C_m", label = "")
[ Plots.plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,10,n], label = "", c = :black) for n in axes(res_p,3) ]

##
p = Plots.plot(plt_CDi_ff, plt_CL_ff, plt_Cm_ff, layout = (1,3), size = (1300, 400))

##
savefig(p, "coeffs.png")