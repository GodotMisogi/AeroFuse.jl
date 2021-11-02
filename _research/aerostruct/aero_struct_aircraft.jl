##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using TimerOutputs
using ComponentArrays

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((2,4,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [0.0, 0.0, 0.0],
            spans     = [5.0, 0.3],
            dihedrals = [0., 45.],
            LE_sweeps = [5., 60.]);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0., 0.08],
             angle     = 3.,
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
                 axis      = [1., 0., 0.]);

## Meshing and assembly

# Wing
wing_n_span   = [8, 3]
wing_n_chord  = 6
vlm_mesh      = chord_coordinates(wing, wing_n_span, wing_n_chord)
cam_mesh      = camber_coordinates(wing, wing_n_span, wing_n_chord)
wing_panels   = make_panels(vlm_mesh)
wing_cambers  = make_panels(cam_mesh)
wing_normals  = panel_normal.(wing_cambers)

# Other panels
htail_panels, htail_normals = panel_wing(htail, 6, 6)
vtail_panels, vtail_normals = panel_wing(vtail, 6, 6)

# Horseshoes
wing_horsies  = Horseshoe.(wing_panels,  wing_normals)
htail_horsies = Horseshoe.(htail_panels,  htail_normals)
vtail_horsies = Horseshoe.(vtail_panels,  vtail_normals)

# Aircraft assembly
aircraft = Dict(
                "Wing"            => wing_horsies,
                "Horizontal Tail" => htail_horsies,
                "Vertical Tail"   => vtail_horsies,
               );

wing_mac = mean_aerodynamic_center(wing);

## Aerodynamic case
ac_name = "My Aircraft"
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
ρ       = 0.98
ref     = [wing_mac[1], 0., 0.]
V, α, β = 25.0, 3.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

## Solve aerodynamic case for initial vector
@time data =
    solve_case(aircraft, fs;
               rho_ref          = ρ,         # Reference density
               r_ref            = ref,       # Reference point for moments
               area_ref         = S,         # Reference area
               span_ref         = b,         # Reference span
               chord_ref        = c,         # Reference chord
               name             = ac_name,   # Aircraft name
               print_components = true,      # Prints the results for each component
              );

## Data collection
Γs = data[ac_name][end]
CFs, CMs, Γ0_wing = data["Wing"][3:end];

## Aerodynamic forces and center locations
vlm_acs    = bound_leg_center.(wing_horsies)
vlm_forces = force.(CFs, dynamic_pressure(ρ, V), S)

# FEM mesh
fem_w    = 0.40
fem_mesh = make_beam_mesh(vlm_mesh, fem_w)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 60 * 9.81
load_factor = 1.3;

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
rs    = range(2e-2, stop = 8e-3, length = length(Ls) ÷ 2)  # Outer radius, m
ts    = range(8e-3, stop = 2e-3, length = length(Ls) ÷ 2)  # Thickness, m
r     = [ reverse(rs); rs ]
t     = [ reverse(ts); ts ]

tubes    = Tube.(Ref(aluminum), Ls, r, t)

## Stiffness matrices, loads and constraints
Ks        = build_big_stiffy(tubes, fem_mesh, vlm_mesh)
cons      = [length(fem_mesh) ÷ 2]
stiffy    = build_stiffness_matrix(Ks, cons)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

## Solve system
dx = solve_cantilever_beam(Ks, fem_loads, cons)
Δx = [ zeros(6); dx[:] ]

## Aerostructural residual
#==========================================================================================#

other_horsies = [ htail_horsies[:]; vtail_horsies[:] ]

# Set up initial guess and function
solve_aerostructural_residual!(R, x) =
    solve_coupled_residual!(R, x,
                            V, deg2rad(β), ρ, Ω, # Aerodynamic state
                            vlm_mesh, cam_mesh,  # Aerodynamic variables
                            other_horsies,       # Other aerodynamic parameters
                            fem_mesh, stiffy,    # Structural variables
                            weight, load_factor) # Load factor variables


# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = Γs,
                    structures   = Δx,
                    load_factor  = deg2rad(α))

## Solve nonlinear system
reset_timer!()
@timeit "Solving Residuals" res_aerostruct =
    nlsolve(solve_aerostructural_residual!, x0,
            method         = :newton,
            show_trace     = true,
            # extended_trace = true,
            autodiff       = :forward,
           );
print_timer()

## Get zero
x_opt = res_aerostruct.zero
Γ_opt = x_opt.aerodynamics
δ_opt = x_opt.structures[7:end]
α_opt = x_opt.load_factor

## Compute displacements
dx  = @views reshape(δ_opt, 6, length(fem_mesh))
dxs = @views SVector.(dx[1,:], dx[2,:], dx[3,:])
Ts  = rotation_matrix(dx[4:6,:])

## Perturb VLM and camber meshes
new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
new_panels   = make_panels(new_vlm_mesh)

new_cam_mesh   = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)
new_cam_panels = make_panels(new_cam_mesh)
new_normals    = panel_normal.(new_cam_panels)

##
Γ_wing      = @views reshape(Γ_opt[1:length(wing_panels)], size(wing_panels))
new_horsies = Horseshoe.(new_panels, new_normals)
all_horsies = [ new_horsies[:]; other_horsies[:] ]

## Aerodynamic forces and center locations
U_opt      = freestream_to_cartesian(-V, α_opt, deg2rad(β))
vlm_acs    = bound_leg_center.(new_horsies)
vlm_forces = nearfield_forces(Γ_wing, new_horsies, Γ_opt, all_horsies, U_opt, Ω, ρ)
all_forces = nearfield_forces(Γ_opt, all_horsies, Γ_opt, all_horsies, U_opt, Ω, ρ)

## New beams and loads
new_fem_mesh = make_beam_mesh(new_vlm_mesh, fem_w)
fem_loads    = compute_loads(vlm_acs, vlm_forces, new_fem_mesh)

## Compute stresses
δxs = eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes
σs  = reduce(hcat, von_mises_stress.(tubes, δxs, δθs))

## Check numbers
lift     = sum(all_forces)[3]
load_fac = lift * cos(α_opt) / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $V m/s")
println("Angle of attack: $(rad2deg(α_opt))ᵒ")

## Generate DataFrame
df = DataFrame(permutedims([ fem_loads; dx ]), :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz, :dx, :dy, :dz, :θx, :θy, :θz])

## Plotting
#==========================================================================================#

# Beam loads and stresses
fem_plot     = reduce(hcat, chop_coordinates(fem_mesh, 1))
new_fem_plot = reduce(hcat, chop_coordinates(new_fem_mesh, 1))
loads_plot   = fem_loads
σs_max       = maximum.(eachcol(σs))
σs_norm      = σs_max ./ maximum(σs_max)
σ_norms      = [ σs_norm; σs_norm[end] ]

# Panels
wing_panel_plot  = plot_panels(wing_panels[:])
htail_panel_plot = plot_panels(htail_panels[:])
vtail_panel_plot = plot_panels(vtail_panels[:])

# Aerodynamic centers and forces
ac_plot    = reduce(hcat, vlm_acs)
force_plot = reduce(hcat, vlm_forces)

# Cambers
cam_plot     = plot_panels(wing_cambers[:])
new_cam_plot = plot_panels(new_cam_panels[:])

# Displacements
new_vlm_mesh_plot = reduce(hcat, new_vlm_mesh)
new_panel_plot    = plot_panels(new_panels[:])

# Axes
xs_plot   = reduce(hcat, (fem_mesh[1:end-1] + fem_mesh[2:end]) / 2)
axes      = axis_transformation(fem_mesh, vlm_mesh)
axes_plot = reshape(reduce(hcat, axes), 3, 3, length(axes))
cs_plot   = axes_plot[:,1,:]
ss_plot   = axes_plot[:,2,:]
ns_plot   = axes_plot[:,3,:]

# Planforms
wing_plan  = plot_wing(wing)
nwing_plan = plot_wing(new_cam_mesh)
htail_plan = plot_wing(htail)
vtail_plan = plot_wing(vtail)

# Streamlines
seed    = chop_coordinates(new_cam_mesh[end,:], 3)
streams = plot_streams(fs, seed, all_horsies, Γ_opt, 5, 100);

## Plot
using Plots
using LaTeXStrings

# gr()
plotlyjs(dpi = 300)
# pyplot(dpi = 300, size = (900, 600))
# pgfplotsx(size = (900, 600))

aircraft_plot =
    plot(xaxis = "x", yaxis = "y", zaxis = "z",
         camera = (-75, 20),
         xlim = (0, b/2),
     #     ylim = (-b/2, b/2),
         zlim = (-b/8, b/4),
         bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :topleft,
         title = "Coupled Aerostructural Analysis"
        )

# Panels
[ plot!(pans, color = :lightgray, label = ifelse(i == 1, "Original Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(cam_plot) ]
[ plot!(pans, color = RGBA(0.5, 0.5, 0.8, 0.5), label = ifelse(i == 1, "Deflected Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_cam_plot) ]
[ plot!(pans, color = :brown, label = :none, linestyle = :solid) for (i, pans) in enumerate(htail_panel_plot) ]
[ plot!(pans, color = :brown, label = :none, linestyle = :solid) for (i, pans) in enumerate(vtail_panel_plot) ]

# Planforms
plot!(wing_plan, color = :gray, label = "Original Wing", linestyle = :solid)
plot!(nwing_plan, color = :blue, label = "Deflected Wing")
plot!(htail_plan, color = :brown, label = "Horizontal Tail")
plot!(vtail_plan, color = :brown, label = "Vertical Tail")

# Beams
thickness = 2.5
r_norm = [ r; r[end]] / maximum(r) * thickness
plot!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_norm)
plot!(new_fem_plot[1,:], new_fem_plot[2,:], new_fem_plot[3,:], m = (thickness, 0.8, :viridis, Plots.stroke(0)), zcolor = σ_norms, cbar = true, label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_norm)

# Streamlines
[ plot!(stream,  color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(streams) ]

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
#         color = :darkcyan, label = :none)
# quiver!(xs_plot[1,:], xs_plot[2,:], xs_plot[3,:],
#         quiver=(ss_plot[1,:], ss_plot[2,:], ss_plot[3,:]) .* 1e-1,
#         color = :black, label = :none)
# quiver!(xs_plot[1,:], xs_plot[2,:], xs_plot[3,:],
#         quiver=(ns_plot[1,:], ns_plot[2,:], ns_plot[3,:]) .* 1e-1,
#         color = :red, label = :none)

# savefig(aircraft_plot, "plots/AerostructAircraft.pdf")
plot!()