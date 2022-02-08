##
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using TimerOutputs

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = WingSection(root_foil  = naca4(2,4,1,2),
                   span       = 2.6,
                   dihedral   = 1.0,
                   sweep        =  20.0,
                   taper      = 0.5,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing),
print_info(wing, "Wing")

# Mesh details
wing_mesh = WingMesh(wing, [12], 6, 
                     span_spacing   = [Cosine()]
                    )

aircraft = ComponentVector(wing = make_horseshoes(wing_mesh));

## Aerodynamic case

# Freestream conditions
fs      = Freestream(alpha = 1.0, 
                     beta  = 0.0, 
                     omega = [0.,0.,0.])

# Reference values
refs    = References(
                     speed    = 15.0,
                     density  = 0.98, 
                     area     = projected_area(wing),   
                     span     = span(wing), 
                     chord    = mean_aerodynamic_chord(wing), 
                     location = [ x_w; 0.; 0. ]
                    )

# Solve aerodynamic case for initial vector
@time data = solve_case(aircraft, fs, refs; 
                        print = true);

## Aerodynamic forces and center locations
vlm_acs    = bound_leg_center.(data.vortices.wing)
vlm_forces = surface_forces(data).wing

## FEM mesh
fem_w    = 0.40
fem_mesh = make_beam_mesh(wing_mesh.chord_mesh, fem_w)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 10 * 9.81
load_factor = 1.0;

## Structural variables

# Material properties
E     = 85e9   # Elastic modulus, N/m²
G     = 25e9   # Shear modulus, N/m²
σ_max = 350e6  # Yield stress with factor of safety 2.5, N/m²
rho   = 1.6e3  # Density, kg/m³
ν     = 0.3    # Poisson ratio (UNUSED FOR NOW)

aluminum = Material(E, G, σ_max, rho)

## Beam properties
Ls    = norm.(diff(fem_mesh))                              # Beam lengths, m
rs    = range(6e-3, stop = 4e-3, length = length(Ls) ÷ 2)  # Outer radius, m
ts    = range(1e-3, stop = 8e-4, length = length(Ls) ÷ 2)  # Thickness, m
r     = [ reverse(rs); rs ]
t     = [ reverse(ts); ts ]

tubes = Tube.(Ref(aluminum), Ls, r, t)

## Stiffness matrix, loads and constraints
Ks        = build_big_stiffy(tubes, fem_mesh, wing_mesh.chord_mesh)
cons      = [length(fem_mesh) ÷ 2]
stiffy    = build_stiffness_matrix(Ks, cons)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

## Solve system
dx = solve_cantilever_beam(Ks, fem_loads, cons)
Δx = ComponentArray(
                    constraint   = zeros(6),
                    displacement = dx
                   )

## Aerostructural residual
#==========================================================================================#

# Set up initial guess and function

function aerostructural_problem(V, β, ρ, Ω, wing_mesh, fem_mesh, stiffy, weight, load_factor)
    f!(R, x) =
        solve_coupled_residual!(R, x,
                                V, β, ρ, Ω,          # Aerodynamic state
                                wing_mesh.chord_mesh,  # VLM mesh
                                wing_mesh.camber_mesh,  # Camber mesh
                                fem_mesh, stiffy,    # Structural variables
                                weight, load_factor) # Load factor variables
end

solve_aerostruct! = aerostructural_problem(refs.speed, fs.beta, refs.density, fs.omega, wing_mesh, fem_mesh, stiffy, weight, load_factor)

# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = data.circulations.wing,
                    structures   = Δx,
                    load_factor  = fs.alpha)

##
# using ForwardDiff
# using Zygote

# function newton_raphson(f!, x0; max_iters = 50, tol = 1e-9)
#     x = copy(x0)
#     R = similar(x)
#     ∂R∂x = Matrix{eltype(x)}(undef, length(R), length(x))
#     ε = 1e5
#     i = 0
#     for i = 1:max_iters
#         ForwardDiff.jacobian!(∂R∂x, f!, R, x)
#         dx   = ∂R∂x \ -R
#         if ε <= tol return x end # Needs NAN checks and everything like NLsolve
#         ε    = LinearAlgebra.norm(dx)
#         @show (i, ε)
#         x  .+= dx
#         i   += 1
#     end
#     return x
# end

# x = @time newton_raphson(solve_aerostruct!, x0)

# x = @time aerostruct_gauss_seidel(x0, V, deg2rad(β), ρ, Ω, wing_mesh.chord_mesh, wing_mesh.camber_mesh, fem_mesh, stiffy, weight, load_factor; max_iters = 50, tol = 1e-9)

##
reset_timer!()
@timeit "Solving Residuals" res_aerostruct =
    nlsolve(solve_aerostruct!, x0,
            method         = :newton,
            autodiff       = :forward,
            # show_trace     = true,
            # extended_trace = true,
           );
print_timer()

## Get zero
x_opt = res_aerostruct.zero
Γ_opt = x_opt.aerodynamics
δ_opt = x_opt.structures.displacement
α_opt = x_opt.load_factor

## Compute displacements
dx  = δ_opt
dxs = mesh_translation(dx)
Ts  = mesh_rotation(dx)

## Perturb VLM and camber meshes
new_chord_mesh = transfer_displacements(dxs, Ts, wing_mesh.chord_mesh, fem_mesh)
new_camber_mesh = transfer_displacements(dxs, Ts, wing_mesh.camber_mesh, fem_mesh)

# new_panels     = make_panels(new_vlm_meh)
# new_cam_panels = make_panels(new_camber_mesh)
# new_normals    = panel_normal.(new_cam_panels)
using Setfield
new_aircraft = ComponentArray(wing = Horseshoe.(make_panels(new_chord_mesh), panel_normal.(make_panels(new_camber_mesh))))
new_fs       = @set fs.alpha = α_opt
system       = solve_case(new_aircraft, new_fs, refs);

## Aerodynamic forces and center locations
U_opt       = body_frame_velocity(new_fs)
vlm_acs     = bound_leg_center.(system.vortices.wing)
vlm_forces  = surface_forces(system).wing

## New beams and loads
new_fem_mesh = make_beam_mesh(new_chord_mesh, fem_w)
fem_loads    = compute_loads(vlm_acs, vlm_forces, new_fem_mesh)

## Compute stresses
δxs = eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes
σs  = combinedimsview(von_mises_stress.(tubes, δxs, δθs))

## Check numbers
 lift     = sum(vlm_forces)[3]
load_fac = lift * cos(α_opt) / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(refs.speed) m/s")
println("Angle of attack: $(rad2deg(α_opt))ᵒ")

## Generate DataFrame
df = DataFrame(permutedims([ fem_loads; dx ]), :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz, :dx, :dy, :dz, :θx, :θy, :θz])

## Fuel burn computation
#==========================================================================================#

## Get coefficients
CD, CY, CL, Cl, Cm, Cn = nearfield(system)

## Compute fuel burn using Bréguet equation
fuel_fractions(ffs) = 1 - prod(ffs)
breguet_fuel_fraction(R, V, SFC, CD, CL) = exp(-R * SFC / (V * (CL / CD)))
breguet_fuel_burn(W, args...) = W * (1 - breguet_fuel_fraction(args...))

SFC = 0.001
R  = 50000
Ws = 0
V  = refs.speed

m_f = breguet_fuel_burn(weight, R, V, SFC, CD, CL)

## Plotting
#==========================================================================================#

# Beam loads and stresses
fem_plot     = combinedimsview(chop_coordinates(fem_mesh, 1))
new_fem_plot = combinedimsview(chop_coordinates(new_fem_mesh, 1))
loads_plot   = fem_loads
σs_max       = maximum.(eachcol(σs))
σs_norm      = σs_max ./ maximum(σs_max)
σ_norms      = [ σs_norm; σs_norm[end] ]

# Aerodynamic centers and forces
cam_panels = make_panels(new_camber_mesh)
ac_plot    = combinedimsview(vlm_acs[:])
force_plot = combinedimsview(vlm_forces[:])

# Cambers
# cam_plot     = plot_panels(cam_panels)
# new_cam_plot = plot_panels(new_cam_panels)

# Displacements
new_chord_mesh_plot = combinedimsview(new_chord_mesh)
# new_panel_plot    = plot_panels(new_panels)

# Axes
xs_plot   = combinedimsview((fem_mesh[1:end-1] + fem_mesh[2:end]) / 2)
# axes      = axis_transformation(fem_mesh, new_chord_mesh)
# axes_plot = combinedimsview(axes)
# cs_plot   = axes_plot[:,1,:]
# ss_plot   = axes_plot[:,2,:]
# ns_plot   = axes_plot[:,3,:]

# Planforms
wing_plan  = plot_planform(wing)
nwing_plan = plot_planform(new_camber_mesh)

# Streamlines
seed    = chop_coordinates(new_camber_mesh[1,:], 3)
streams = streamlines(new_fs, refs, seed, system.vortices, Γ_opt, 1, 100);

## Plot
using LaTeXStrings
# using GeometryTypes

# using CairoMakie
# CairoMakie.activate!()
using GLMakie
# activate!()

const LS = LaTeXString

## Surface pressure coefficients
CFs, CMs         = surface_coefficients(data; axes = Wind())
new_CFs, new_CMs = surface_coefficients(system; axes = Wind())

cps     = norm.(CFs) * S
new_cps = norm.(new_CFs) * S

wing_cp_points     = extrapolate_point_mesh(cps.wing)
new_wing_cp_points = extrapolate_point_mesh(new_cps.wing)

## Spanwise loads
cl_loading(Γs, V, c) = vec(sum(Γs, dims = 1)) / (0.5 * V * c)

hs_pts      = horseshoe_point.(data.vortices)
wing_ys     = @view combinedimsview(hs_pts.wing[1,:])[2,:]

new_hs_pts  = horseshoe_point.(system.vortices)
new_wing_ys = @view combinedimsview(new_hs_pts.wing[1,:])[2,:]

wing_ll     = spanwise_loading(chord_panels(wing_mesh), CFs.wing, S)

new_wing_ll = spanwise_loading(make_panels(new_chord_mesh), new_CFs.wing, S)

## Mesh connectivities
triangle_connectivities(inds) = @views [ inds[1:end-1,1:end-1][:] inds[1:end-1,2:end][:]   inds[2:end,2:end][:]   ;
                                           inds[2:end,2:end][:]   inds[2:end,1:end-1][:] inds[1:end-1,1:end-1][:] ]

wing_cam_connec     = triangle_connectivities(LinearIndices(wing_mesh.camber_mesh))
new_wing_cam_connec = triangle_connectivities(LinearIndices(new_camber_mesh))

## Figure plot
fig   = Figure(resolution = (1280, 720))

scene = LScene(fig[1:4,1])
ax    = fig[1:4,2] = GridLayout()

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
plot_spanload!(ax, new_wing_ll, LS("Deflected Wing"))

# Meshes
m1 = poly!(scene, wing_mesh.camber_mesh[:], wing_cam_connec, color = wing_cp_points[:], shading = true, transparency = true)
m2 = poly!(scene, new_camber_mesh[:], new_wing_cam_connec, color = new_wing_cp_points[:])

# Borders
# lines!(scene, plot_planform(wing))
l1 = [ lines!(scene, pts, color = :grey, transparency = true) for pts in plot_panels(make_panels(wing_mesh.camber_mesh)) ]
l2 = [ lines!(scene, pts, color = :grey, transparency = true) for pts in plot_panels(make_panels(new_camber_mesh)) ]

# Streamlines
[ lines!(scene, stream[:], color = :lightblue) for stream in eachcol(streams) ]

cam = cam3d!(scene.scene)

leg = Legend(ax[4,1:2], ax_cl)
leg.tellwidth = true

fig[0, :] = Label(fig, LS("Aerostructural Analysis"), textsize = 20)

fig
##
save("plots/AerostructWingMakie.pdf", fig)

## Plot
# using Plots
# using LaTeXStrings

# # pyplot(dpi = 300)
# # pgfplotsx(size = (900, 600))

# aircraft_plot =
#     plot(xaxis = L"$x$", yaxis = L"$y$", zaxis = L"$z$",
#          camera = (-85, 20),
#          xlim = (-b/4, 3b/4),
#      #     ylim = (-b/2, b/2),
#          zlim = (-b/8, b/4),
#          bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
#          legend = :bottomright,
#          title = "Coupled Aerostructural Analysis"
#         )

# # Panels
# [ plot!(pans, color = :lightgray, label = ifelse(i == 1, "Original Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(cam_plot) ]
# [ plot!(pans, color = RGBA(0.5, 0.5, 0.8, 0.7), label = ifelse(i == 1, "Deflected Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_cam_plot) ]

# # Planform
# plot!(wing_plan,  color = :gray, label = "Original Wing Planform", linestyle = :solid)
# plot!(nwing_plan, color = :blue, label = "Deflected Wing Planform")

# # Beams
# thickness = 2.5
# r_norm = [ r; r[end]] / maximum(r) * thickness
# plot!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_norm)
# plot!(new_fem_plot[1,:], new_fem_plot[2,:], new_fem_plot[3,:], color = RGBA.(σ_norms, 0.5, 0.6, 1.0), label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_norm)

# # Streamlines
# [ plot!(stream,  color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(eachcol(streams)) ]

# Forces
# quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
#         quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1,
#         label = "Panel Forces", color = :orange)
# scatter!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], label = "Aerodynamic Centers")
# quiver!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:],
#         quiver=(loads_plot[1,:], loads_plot[2,:], loads_plot[3,:] ) .* 0.1,
#         label = "Beam Forces")

# Axis systems
# quiver!(xs_plot[1,:], xs_plot[2,:], xs_plot[3,:],
#         quiver=(cs_plot[1,:], cs_plot[2,:], cs_plot[3,:]) .* 1e-1,
#         color = :orange, label = :none)
# quiver!(xs_plot[1,:], xs_plot[2,:], xs_plot[3,:],
#         quiver=(ss_plot[1,:], ss_plot[2,:], ss_plot[3,:]) .* 1e-1,
#         color = :black, label = :none)
# quiver!(xs_plot[1,:], xs_plot[2,:], xs_plot[3,:],
#         quiver=(ns_plot[1,:], ns_plot[2,:], ns_plot[3,:]) .* 1e-1,
#         color = :red, label = :none)

# savefig(aircraft_plot, "plots/AerostructWing.pdf")
# plot!()