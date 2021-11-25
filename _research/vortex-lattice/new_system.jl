## Wing analysis case
using Revise
using AeroMDAO
using ComponentArrays
using LinearAlgebra

## Surfaces

# Wing
wing = Wing(foils     = Foil.(fill(naca4(2,4,1,2), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            LE_sweeps = [5.]);

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

## Meshing

## Wing
wing_n_span   = [12]
wing_n_chord  = 6
wing_vlm_mesh = chord_coordinates(wing, wing_n_span, wing_n_chord)
wing_cam_mesh = camber_coordinates(wing, wing_n_span, wing_n_chord)
wing_panels   = make_panels(wing_vlm_mesh)
wing_cambers  = make_panels(wing_cam_mesh)
wing_normals  = panel_normal.(wing_cambers)
wing_horsies  = Horseshoe.(wing_panels,  wing_normals)

struct WingMesh{M <: AbstractWing, N <: Integer, P} <: AbstractWing
    surf    :: M
    n_span  :: N
    n_chord :: P
end

WingMesh(surf :: Union{Wing{T}, HalfWing{T}}, n_span :: AbstractVector{M}, n_chord :: N) where {T <: Real, N <: Integer, M <: Integer} = VLMWing{T,N,M}(surf, n_span, n_chord) 

# AeroMDAO.AircraftGeometry.chord_coordinates(wing :: VLMWing) = chord_coordinates(wing.surf, wing.n_span, wing.n_chord)

# vlm_wing = WingMesh(wing, wing_n_span, wing_n_chord)

# Horizontal tail
htail_n_span   = [6]
htail_n_chord  = 6
htail_vlm_mesh = chord_coordinates(htail, htail_n_span, htail_n_chord)
htail_cam_mesh = camber_coordinates(htail, htail_n_span, htail_n_chord)
htail_panels   = make_panels(htail_vlm_mesh)
htail_cambers  = make_panels(htail_cam_mesh)
htail_normals  = panel_normal.(htail_cambers)
htail_horsies  = Horseshoe.(htail_panels, htail_normals)

# Vertical tail
vtail_n_span   = [6]
vtail_n_chord  = 6
vtail_vlm_mesh = chord_coordinates(vtail, vtail_n_span, vtail_n_chord)
vtail_cam_mesh = camber_coordinates(vtail, vtail_n_span, vtail_n_chord)
vtail_panels   = make_panels(vtail_vlm_mesh)
vtail_cambers  = make_panels(vtail_cam_mesh)
vtail_normals  = panel_normal.(vtail_cambers)
vtail_horsies  = Horseshoe.(vtail_panels, vtail_normals)

aircraft_panels = ComponentArray(
                                 wing  = wing_panels,
                                 htail = htail_panels, 
                                 vtail = vtail_panels
                                )

aircraft = ComponentArray(
                          wing  = wing_horsies,
                          htail = htail_horsies,
                          vtail = vtail_horsies
                         )

wing_mac = mean_aerodynamic_center(wing)
x_w      = wing_mac[1]

## Case
ac_name = "My Aircraft"
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
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
                     );

    CFs, CMs = surface_coefficients(data; axes = Wind())
    FFs = farfield_coefficients(data)
end;

## Spanwise forces
function span_loading(panels, CFs, Γs, V, alpha, beta, c)
    CFs  = body_to_wind_axes.(CFs, alpha, beta)
    CDis = @. getindex(CFs, 1)
    CYs  = @. getindex(CFs, 2)
    CLs  = @. getindex(CFs, 3)

    area_scale  = S ./ sum(panel_area, panels, dims = 1)[:]
    span_CDis   = sum(CDis, dims = 1)[:] .* area_scale
    span_CYs    = sum(CYs,  dims = 1)[:] .* area_scale
    span_CLs    = sum(CLs,  dims = 1)[:] .* area_scale
    CL_loadings = sum(Γs,   dims = 1)[:] / (0.5 * V * c)

    span_CDis, span_CYs, span_CLs, CL_loadings
end

span_loading(panels, CFs, Γs, fs, c) = span_loading(panels, CFs, Γs, fs.V, fs.alpha, fs.beta, c)

hs_pts   = horseshoe_point.(data.horseshoes)
wing_ys  = getindex.(hs_pts.wing[1,:], 2)
htail_ys = getindex.(hs_pts.htail[1,:], 2)
vtail_ys = getindex.(hs_pts.vtail[1,:], 2)

wing_CDis, wing_CYs, wing_CLs, wing_CL_loadings = span_loading(wing_panels, CFs.wing, data.circulations.wing, fs, c)
htail_CDis, htail_CYs, htail_CLs, htail_CL_loadings = span_loading(htail_panels, CFs.htail, data.circulations.htail, fs, c)
vtail_CDis, vtail_CYs, vtail_CLs, vtail_CL_loadings = span_loading(vtail_panels, CFs.vtail, data.circulations.vtail, fs, c);

## Plotting
using GLMakie
using LaTeXStrings

const LS = LaTeXString

## Streamlines
# Spanwise distribution
span_points = 20
init        = chop_leading_edge(wing, span_points)
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy,  dz])
                init .+ Ref([dx, dy, -dz]) ];

distance = 5
num_stream_points = 100
streams = plot_streams(fs, seed, data.horseshoes, data.circulations, distance, num_stream_points);

## Mesh connectivities
triangle_connectivities(inds) = @views [ inds[1:end-1,1:end-1][:] inds[1:end-1,2:end][:]   inds[2:end,2:end][:]   ;
                                           inds[2:end,2:end][:]   inds[2:end,1:end-1][:] inds[1:end-1,1:end-1][:] ]

wing_cam_connec  = triangle_connectivities(LinearIndices(wing_cam_mesh))
htail_cam_connec = triangle_connectivities(LinearIndices(htail_cam_mesh))
vtail_cam_connec = triangle_connectivities(LinearIndices(vtail_cam_mesh));

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

## Surfave velocities
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

## Coordinates
fig  = Figure()

scene = LScene(fig[1:4,1])
ax1   = fig[1,2] = GLMakie.Axis(fig, ylabel = L"C_{D_i}", title = LS("Spanload Distributions"))
ax2   = fig[2,2] = GLMakie.Axis(fig, ylabel = L"C_Y",)
ax3   = fig[3,2] = GLMakie.Axis(fig, xlabel = L"y", ylabel = L"C_L")

# Meshes
m1 = poly!(scene, wing_cam_mesh[:],  wing_cam_connec,  color =  wing_cp_points[:])
m2 = poly!(scene, htail_cam_mesh[:], htail_cam_connec, color = htail_cp_points[:])
m3 = poly!(scene, vtail_cam_mesh[:], vtail_cam_connec, color = vtail_cp_points[:])
l1 = lines!.(scene, streams, color = :green)

# wing_surf = surface_coordinates(wing, [30], 60)
# surf_connec = square_connectivities(LinearIndices(wing_surf))
# wing_surf_mesh = mesh(wing_surf[:], surf_connec)
# w1 = wireframe!(scene, wing_surf_mesh.plot[1][], color = :grey, alpha = 0.1)

# Spanload plot
function plot_spanload!(fig, ys, CDis, CYs, CLs, CL_loadings, name = "Wing")
    lines!(fig[1,2], ys, CDis, label = name,)
    lines!(fig[2,2], ys, CYs, label = name,)
    lines!(fig[3,2], ys, CLs, label = name,)
    # lines!(fig[3,2], ys, CL_loadings, label = "$name Loading")

    nothing
end

plot_spanload!(fig, wing_ys, wing_CDis, wing_CYs, wing_CLs, wing_CL_loadings, LS("Wing"))
plot_spanload!(fig, htail_ys, htail_CDis, htail_CYs, htail_CLs, htail_CL_loadings, LS("Horizontal Tail"))
plot_spanload!(fig, vtail_ys, vtail_CDis, vtail_CYs, vtail_CLs, vtail_CL_loadings, LS("Vertical Tail"))

# Legend
axl = fig[4,2] = GridLayout()
Legend(fig[4,1:2], ax3)
fig[0, :] = Label(fig, LS("Vortex Lattice Analyses"))

fig.scene

##
# hs_pts = Tuple.(bound_leg_center.(horses))[:]
# arrows!(scene, getindex.(hs_pts, 1), getindex.(hs_pts, 2), getindex.(hs_pts, 3), 
#                 CDis[:], CYs[:], CLs[:], 
#                 arrowsize = Vec3f.(0.3, 0.3, 0.4),
#                 lengthscale = 10,
#                 label = "Forces (Exaggerated)")