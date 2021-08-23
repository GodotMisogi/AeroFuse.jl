##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using TimerOutputs
using ComponentArrays
using SparseArrays

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((2,4,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [0.0, 0.0, 0.0],
            spans     = [5.0, 0.3],
            dihedrals = [0., 45.],
            sweep_LEs = [5., 60.]);

# Horizontal tail 
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39],
             position  = [4., 0., 0.],
             angle     = 0.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)), 
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.]);

## Meshing and assembly
span_wing, chord_wing       = [8, 3], 6
wing_panels,  wing_normals  = panel_wing(wing, span_wing, chord_wing;)

span_htail, chord_htail     = [6], 6
htail_panels, htail_normals = panel_wing(htail, span_htail, chord_htail;)

vtail_panels, vtail_normals = panel_wing(vtail, 6, 6;)

# Aircraft assembly
aircraft = Dict(
                "Wing"            => (wing_panels,  wing_normals),
                "Horizontal Tail" => (htail_panels, htail_normals),
                "Vertical Tail"   => (vtail_panels, vtail_normals),
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
               rho_ref     = ρ, 		# Reference density
               r_ref       = ref, 		# Reference point for moments
               area_ref    = S, 		# Reference area
               span_ref    = b, 		# Reference span
               chord_ref   = c, 		# Reference chord
               name        = ac_name,	# Aircraft name
               print_components = true,	# Prints the results for each component
              );

## Data collection
Γs = data[ac_name][end]
CFs_wing, CMs_wing, hs_wing, Γ0_wing = data["Wing"][3:end];
CFs_htail, CMs_htail, hs_htail, Γ0_htail = data["Horizontal Tail"][3:end];
Γ0_vtail = data["Vertical Tail"][end];

## Wing FEM setup
vlm_mesh_wing   = chord_coordinates(wing, span_wing, chord_wing)
cam_mesh_wing   = camber_coordinates(wing, span_wing, chord_wing)

cam_panels_wing = make_panels(cam_mesh_wing)

vlm_acs_wing    = bound_leg_center.(hs_wing)
vlm_forces_wing = force.(CFs_wing, dynamic_pressure(ρ, V), S)

fem_weight_wing = 0.40
fem_mesh_wing   = make_beam_mesh(vlm_mesh_wing, fem_weight_wing)

aluminum = Material(       # Aluminum properties
                    85e9,  # Elastic modulus, N/m² 
                    25e9,  # Shear modulus, N/m², 
                    350e6, # Yield stress with factor of safety 2.5, N/m², 
                    1.6e3, # Density, kg/m³
                    )

Ls_wing = norm.(diff(fem_mesh_wing))                              # Beam lengths, m 
rs_wing = range(2e-2, stop = 8e-3, length = length(Ls_wing) ÷ 2)  # Outer radius, m
ts_wing = range(8e-3, stop = 2e-3, length = length(Ls_wing) ÷ 2)  # Thickness, m
r_wing  = [ reverse(rs_wing); rs_wing ]
t_wing  = [ reverse(ts_wing); ts_wing ]

tubes_wing     = Tube.(Ref(aluminum), Ls_wing, r_wing, t_wing)
Ks_wing        = build_big_stiffy(tubes_wing, fem_mesh_wing, vlm_mesh_wing)
cons_wing      = [length(fem_mesh_wing) ÷ 2]
stiffy_wing    = build_stiffness_matrix(Ks_wing, cons_wing)
fem_loads_wing = compute_loads(vlm_acs_wing, vlm_forces_wing, fem_mesh_wing)

dx_wing = solve_cantilever_beam(Ks_wing, fem_loads_wing, cons_wing)
Δx_wing = [ zeros(6); dx_wing[:] ]

## Horizontal tail FEM setup
vlm_mesh_htail   = chord_coordinates(htail, span_htail, chord_htail)
cam_mesh_htail   = camber_coordinates(htail, span_htail, chord_htail)

cam_panels_htail = make_panels(cam_mesh_htail)

vlm_acs_htail    = bound_leg_center.(hs_htail)
vlm_forces_htail = force.(CFs_htail, dynamic_pressure(ρ, V), S)

fem_weight_htail = 0.40
fem_mesh_htail   = make_beam_mesh(vlm_mesh_htail, fem_weight_htail)

# Beam properties
Ls_htail = norm.(diff(fem_mesh_htail))                              # Beam lengths, m 
rs_htail = range(8e-3, stop = 2e-3, length = length(Ls_htail) ÷ 2)  # Outer radius, m
ts_htail = range(4e-3, stop = 1e-3, length = length(Ls_htail) ÷ 2)  # Thickness, m
r_htail  = [ reverse(rs_htail); rs_htail ]
t_htail  = [ reverse(ts_htail); ts_htail ]

tubes_htail     = Tube.(Ref(aluminum), Ls_htail, r_htail, t_htail)
Ks_htail        = build_big_stiffy(tubes_htail, fem_mesh_htail, vlm_mesh_htail)
cons_htail      = [length(fem_mesh_htail) ÷ 2]
stiffy_htail    = build_stiffness_matrix(Ks_htail, cons_htail)
fem_loads_htail = compute_loads(vlm_acs_htail, vlm_forces_htail, fem_mesh_htail)

dx_htail = solve_cantilever_beam(Ks_htail, fem_loads_htail, cons_htail)
Δx_htail = [ zeros(6); dx_htail[:] ]

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 60 * 9.81
load_factor = 1.3;

##
stiffy = blockdiag(stiffy_wing, stiffy_htail)

## Aerostructural residual
#==========================================================================================#

other_horsies = Horseshoe.(vtail_panels[:], vtail_normals[:])

vlm_meshes  = [ vlm_mesh_wing, vlm_mesh_htail ]
cam_meshes  = [ cam_mesh_wing, cam_mesh_htail ]
fem_meshes  = [ fem_mesh_wing, fem_mesh_htail ]
fem_weights = [ fem_weight_wing, fem_weight_htail ]
syms        = [ :wing, :htail ]

# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = (wing = Γ0_wing, htail = Γ0_htail, vtail = Γ0_vtail),
                    structures   = (wing = Δx_wing, htail = Δx_htail),
                    load_factor  = deg2rad(α))

# Set up initial guess and function
solve_aerostructural_residual!(R, x) = 
    solve_coupled_residual!(R, x,
                            V, deg2rad(β), ρ, Ω,
                            syms, vlm_meshes, cam_meshes, fem_meshes, 
                            other_horsies, stiffy, weight, load_factor)

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
δ_opt = x_opt.structures
α_opt = x_opt.load_factor

## Compute displacements
Δs    = map((key, n) -> reshape(δ_opt[key][7:end], 6, n), valkeys(δ_opt), length.(fem_meshes))
dx_Ts = translations_and_rotations.(Δs)
dxs   = getindex.(dx_Ts, 1)
Ts    = getindex.(dx_Ts, 2)

## New VLM variables
new_vlm_meshes = transfer_displacements.(dxs, Ts, vlm_meshes, fem_meshes)
new_panels     = make_panels.(new_vlm_meshes)

new_cam_meshes = transfer_displacements.(dxs, Ts, cam_meshes, fem_meshes)
new_cam_panels = make_panels.(new_cam_meshes) 

##
new_horsies = new_horseshoes.(dxs, Ts, vlm_meshes, cam_meshes, fem_meshes)
all_horsies = [ reduce(vcat, vec.(new_horsies)); other_horsies ]

## Aerodynamic forces and center locations
U_opt      = freestream_to_cartesian(-V, α_opt, deg2rad(β))
new_acs    = map(horsies -> bound_leg_center.(horsies), new_horsies)
all_forces = nearfield_forces(Γ_opt, all_horsies, Γ_opt, all_horsies, U_opt, Ω, ρ)

new_Γs     = getindex.(Ref(Γ_opt), syms)
new_forces = nearfield_forces.(new_Γs, new_horsies, Ref(Γs), Ref(all_horsies), Ref(U_opt), Ref(Ω), Ref(ρ))

## New beams and loads
new_fem_meshes = make_beam_mesh.(new_vlm_meshes, fem_weights)
fem_loads      = compute_loads.(new_acs, new_forces, new_fem_meshes)

## Compute stresses
# δxs = eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
# δθs = eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes 
# σs  = reduce(hcat, von_mises_stress.(tubes, δxs, δθs))

## Check numbers
lift     = sum(all_forces)[3]
load_fac = lift * cos(α_opt) / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $V m/s")
println("Angle of attack: $(rad2deg(α_opt))ᵒ")

## Generate DataFrame
dfs = DataFrame.((permutedims ∘ vcat).(fem_loads, Δs), :auto)
rename!.(dfs, Ref([:Fx, :Fy, :Fz, :Mx, :My, :Mz, :dx, :dy, :dz, :θx, :θy, :θz]))

## Plotting
#==========================================================================================#

# Beam loads and stresses
fem_plot     = @. reduce(hcat, chop_coordinates(fem_meshes, 1))
new_fem_plot = @. reduce(hcat, chop_coordinates(new_fem_meshes, 1))
# loads_plot   = fem_loads

# Panels
wing_panel_plot  = plot_panels(wing_panels[:])
htail_panel_plot = plot_panels(htail_panels[:]) 
vtail_panel_plot = plot_panels(vtail_panels[:])

# Aerodynamic centers and forces
ac_plot    = @. reduce(hcat, new_acs)
force_plot = @. reduce(hcat, new_forces)

# Cambers
cam_panels   = @. make_panels(cam_meshes)
cam_plot     = plot_panels(reduce(vcat, vec.(cam_panels)))
new_cam_plot = plot_panels(reduce(vcat, vec.(new_cam_panels)))

# Displacements
new_vlm_mesh_plot = @. reduce(hcat, new_vlm_meshes)
new_panel_plot    = plot_panels(reduce(vcat, vec.(make_panels.(new_vlm_meshes))))

# Planforms
wing_plan   = plot_wing(wing)
htail_plan  = plot_wing(htail)                      
vtail_plan  = plot_wing(vtail)

# New planforms
nwing_plan  = plot_wing(new_cam_meshes[1])
nhtail_plan = plot_wing(new_cam_meshes[2])

# Streamlines
seed    = chop_coordinates(new_cam_meshes[1][end,:], 2)
streams = plot_streams(fs, seed, all_horsies, Γ_opt, 5, 100);

## Plot
using Plots
using LaTeXStrings

pyplot(dpi = 300)
# pgfplotsx(size = (900, 600))

aircraft_plot = 
    plot(xaxis = L"$x$", yaxis = L"$y$", zaxis = L"$z$",
         camera = (-60, 30), 
         xlim = (-b/4, 3b/4),
     #     ylim = (-b/2, b/2),
         zlim = (-b/8, b/8),
         bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :bottomright,
         title = "Coupled Aerostructural Analysis"
        )

# Panels
[ plot!(pans, color = :lightgray, label = ifelse(i == 1, "Original Panels", :none), linestyle = :solid) for (i, pans) in enumerate(cam_plot) ]
[ plot!(pans, color = RGBA(0.5, 0.5, 0.8, 0.5), label = ifelse(i == 1, "Deflected Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_cam_plot) ]
[ plot!(pans, color = :brown, label = :none, linestyle = :solid) for (i, pans) in enumerate(vtail_panel_plot) ]

# Planforms
plot!(wing_plan, color = :gray, label = "Original Wing", linestyle = :solid)
plot!(nwing_plan, color = :blue, label = "Deflected Wing")
plot!(htail_plan, color = :gray, label = "Horizontal Tail")
plot!(nhtail_plan, color = :blue, label = "Deflected Horizontal Tail")
plot!(vtail_plan, color = :brown, label = "Vertical Tail")

# Beams
# thickness = 2.5
# r_norm = [ r; r[end]] / maximum(r) * thickness
# plot!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_norm)
# plot!(new_fem_plot[1,:], new_fem_plot[2,:], new_fem_plot[3,:], color = RGBA.(σ_norms, 0.5, 0.6, 1.0), label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_norm)

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

# savefig(aircraft_plot, "plots/AerostructWingTail.pdf")
plot!()