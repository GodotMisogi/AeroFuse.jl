## TOTALLY not a ripoff of MIT's Dawn Solar HALE aircraft
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using TimerOutputs
using ComponentArrays
using SparseArrays

# Aerostructural analysis case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((2,4,1,2)), 3)),
            chords    = [1.0, 1.0, 0.6],
            twists    = [0.0, 0.0, 0.0],
            spans     = [4.0, 3.0],
            dihedrals = [0., 0.],
            LE_sweeps = [0., 5.]);

print_info(wing, "Lawn Polar Wing")

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0., 0.],
             angle     = 0.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail_u = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                   chords    = [0.7, 0.25],
                   twists    = [0.0, 0.0],
                   spans     = [1.0],
                   dihedrals = [0.],
                   LE_sweeps = [7.97],
                   position  = [4.7, 0, 0],
                   angle     = 90.,
                   axis      = [1., 0., 0.]);

vtail_d = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                   chords    = [0.7, 0.42],
                   twists    = [0.0, 0.0],
                   spans     = [0.4],
                   dihedrals = [0.],
                   LE_sweeps = [7.97],
                   position  = [4.7, 0, 0],
                   angle     = 90.,
                   axis      = [1., 0., 0.]);

vtail = Wing(vtail_d, vtail_u)

# Tailerons
atail_l = Wing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
               chords    = [0.4, 0.2],
               twists    = [0.0, 0.0],
               spans     = [0.5],
               dihedrals = [0.],
               LE_sweeps = [8.],
               position  = [2., 2.5, 0.],
               angle     = 0.,
               axis      = [0., 1., 0.])

atail_r = Wing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
               chords    = [0.4, 0.2],
               twists    = [0.0, 0.0],
               spans     = [0.5],
               dihedrals = [0.],
               LE_sweeps = [8.],
               position  = [2., -2.5, 0.],
               angle     = 0.,
               axis      = [0., 1., 0.]);

## Meshing and assembly

# Wing
wing_n_span    = [6, 6]
wing_n_chord   = 6
wing_vlm_mesh  = chord_coordinates(wing, wing_n_span, wing_n_chord;
                                    # spacings = "cosine"
                                  )
wing_cam_mesh  = camber_coordinates(wing, wing_n_span, wing_n_chord;
                                    # spacings = "cosine"
                                   )
wing_panels    = make_panels(wing_vlm_mesh)
wing_cambers   = make_panels(wing_cam_mesh)
wing_normals   = panel_normal.(wing_cambers)
wing_horsies   = Horseshoe.(wing_panels, wing_normals)

# Horizontal tail
htail_n_span   = [6]
htail_n_chord  = 4
htail_vlm_mesh = chord_coordinates(htail, htail_n_span, htail_n_chord)
htail_cam_mesh = camber_coordinates(htail, htail_n_span, htail_n_chord)
htail_panels   = make_panels(htail_vlm_mesh)
htail_cambers  = make_panels(htail_cam_mesh)
htail_normals  = panel_normal.(htail_cambers)
htail_horsies  = Horseshoe.(htail_panels, htail_normals)

# Vertical tail
vtail_n_span   = [6]
vtail_n_chord  = 4
vtail_vlm_mesh = chord_coordinates(vtail, vtail_n_span, vtail_n_chord)
vtail_cam_mesh = camber_coordinates(vtail, vtail_n_span, vtail_n_chord)
vtail_panels   = make_panels(vtail_vlm_mesh)
vtail_cambers  = make_panels(vtail_cam_mesh)
vtail_normals  = panel_normal.(vtail_cambers)
vtail_horsies  = Horseshoe.(vtail_panels, vtail_normals)

# Tailerons
atail_n_span     = [4]
atail_n_chord    = 2

atail_l_vlm_mesh = chord_coordinates(atail_l, atail_n_span, atail_n_chord)
atail_l_cam_mesh = camber_coordinates(atail_l, atail_n_span, atail_n_chord)
atail_l_panels   = make_panels(atail_l_vlm_mesh)
atail_l_cambers  = make_panels(atail_l_cam_mesh)
atail_l_normals  = panel_normal.(atail_l_cambers)
atail_l_horsies  = Horseshoe.(atail_l_panels, atail_l_normals)

atail_r_vlm_mesh = chord_coordinates(atail_r, atail_n_span, atail_n_chord)
atail_r_cam_mesh = camber_coordinates(atail_r, atail_n_span, atail_n_chord)
atail_r_panels   = make_panels(atail_r_vlm_mesh)
atail_r_cambers  = make_panels(atail_r_cam_mesh)
atail_r_normals  = panel_normal.(atail_r_cambers)
atail_r_horsies  = Horseshoe.(atail_r_panels, atail_r_normals)

# Aircraft assembly
aircraft = ComponentArray(
                          wing    = wing_horsies,
                          htail   = htail_horsies,
                          vtail   = vtail_horsies,
                          atail_l = atail_l_horsies,
                          atail_r = atail_r_horsies,
                         );

## Aerodynamic analsis
#==========================================================================================#

# Reference values
ac_name  = "Lawn Polar"
wing_mac = mean_aerodynamic_center(wing);
S, b, c  = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);
ρ        = 0.98
ref      = [wing_mac[1], 0., 0.]
V, α, β  = 25.0, 0.0, 0.0
Ω        = zeros(3)
fs       = Freestream(V, α, β, Ω)
q        = dynamic_pressure(ρ, V)

## Solve aerodynamic case for initial vector
@time data =
    solve_case(aircraft, fs;
               rho_ref          = ρ,        # Reference density
               r_ref            = ref,      # Reference point for moments
               area_ref         = S,        # Reference area
               span_ref         = b,        # Reference span
               chord_ref        = c,        # Reference chord
               name             = ac_name,  # Aircraft name
               print            = true,     # Print results
               print_components = true,     # Prints the results for each component
              );

## Wing FEM setup
vlm_acs_wing    = bound_leg_center.(data.wing.horseshoes)
vlm_forces_wing = force.(data.wing.CFs, q, S)

wing_beam_ratio = 0.40
wing_fem_mesh   = make_beam_mesh(wing_vlm_mesh, wing_beam_ratio)

aluminum = Material(       # Aluminum properties
                    85e9,  # Elastic modulus, N/m²
                    25e9,  # Shear modulus, N/m²,
                    350e6, # Yield stress with factor of safety 2.5, N/m²,
                    1.6e3, # Density, kg/m³
                   )

Ls_wing = norm.(diff(wing_fem_mesh))                              # Beam lengths, m
rs_wing = range(5e-2, stop = 1e-2, length = length(Ls_wing) ÷ 2)  # Outer radius, m
ts_wing = range(1e-2, stop = 6e-3, length = length(Ls_wing) ÷ 2)  # Thickness, m
r_wing  = [ reverse(rs_wing); rs_wing ]
t_wing  = [ reverse(ts_wing); ts_wing ]

tubes_wing     = Tube.(Ref(aluminum), Ls_wing, r_wing, t_wing)
Ks_wing        = build_big_stiffy(tubes_wing, wing_fem_mesh, wing_vlm_mesh)
cons_wing      = [length(wing_fem_mesh) ÷ 2]
stiffy_wing    = build_stiffness_matrix(Ks_wing, cons_wing)
fem_loads_wing = compute_loads(vlm_acs_wing, vlm_forces_wing, wing_fem_mesh)

dx_wing = solve_cantilever_beam(Ks_wing, fem_loads_wing, cons_wing)
Δx_wing = [ zeros(6); dx_wing[:] ]

## Horizontal tail FEM setup
vlm_acs_htail    = bound_leg_center.(data.htail.horseshoes)
vlm_forces_htail = force.(data.htail.CFs, q, S)

htail_beam_ratio = 0.35
htail_fem_mesh   = make_beam_mesh(htail_vlm_mesh, htail_beam_ratio)

# Beam properties
Ls_htail = norm.(diff(htail_fem_mesh))                              # Beam lengths, m
rs_htail = range(8e-3, stop = 2e-3, length = length(Ls_htail) ÷ 2)  # Outer radius, m
ts_htail = range(6e-4, stop = 2e-4, length = length(Ls_htail) ÷ 2)  # Thickness, m
r_htail  = [ reverse(rs_htail); rs_htail ]
t_htail  = [ reverse(ts_htail); ts_htail ]

tubes_htail     = Tube.(Ref(aluminum), Ls_htail, r_htail, t_htail)
Ks_htail        = build_big_stiffy(tubes_htail, htail_fem_mesh, htail_vlm_mesh)
cons_htail      = [length(htail_fem_mesh) ÷ 2]
stiffy_htail    = build_stiffness_matrix(Ks_htail, cons_htail)
fem_loads_htail = compute_loads(vlm_acs_htail, vlm_forces_htail, htail_fem_mesh)

dx_htail = solve_cantilever_beam(Ks_htail, fem_loads_htail, cons_htail)
Δx_htail = [ zeros(6); dx_htail[:] ]

## Vertical tail FEM setup
vlm_acs_vtail    = bound_leg_center.(data.vtail.horseshoes)
vlm_forces_vtail = force.(data.vtail.CFs, q, S)

vtail_beam_ratio = 0.35
vtail_fem_mesh   = make_beam_mesh(vtail_vlm_mesh, vtail_beam_ratio)

# Beam properties
Ls_vtail = norm.(diff(vtail_fem_mesh))                          # Beam lengths, m
rs_vtail = range(8e-3, stop = 2e-3, length = length(Ls_vtail))  # Outer radius, m
ts_vtail = range(6e-4, stop = 2e-4, length = length(Ls_vtail))  # Thickness, m
r_vtail  = rs_vtail
t_vtail  = ts_vtail

tubes_vtail     = Tube.(Ref(aluminum), Ls_vtail, r_vtail, t_vtail)
Ks_vtail        = build_big_stiffy(tubes_vtail, vtail_fem_mesh, vtail_vlm_mesh)
cons_vtail      = [1]
stiffy_vtail    = build_stiffness_matrix(Ks_vtail, cons_vtail)
fem_loads_vtail = compute_loads(vlm_acs_vtail, vlm_forces_vtail, vtail_fem_mesh)

dx_vtail = solve_cantilever_beam(Ks_vtail, fem_loads_vtail, cons_vtail)
Δx_vtail = [ zeros(6); dx_vtail[:] ]

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 65 * 9.81
load_factor = 1.3;

##
stiffy = blockdiag(stiffy_wing, stiffy_htail, stiffy_vtail)

## Aerostructural residual
#==========================================================================================#

vlm_meshes  = [ wing_vlm_mesh, htail_vlm_mesh, vtail_vlm_mesh ]
cam_meshes  = [ wing_cam_mesh, htail_cam_mesh, vtail_cam_mesh ]
fem_meshes  = [ wing_fem_mesh, htail_fem_mesh, vtail_fem_mesh ]
fem_weights = [ wing_beam_ratio, htail_beam_ratio, vtail_beam_ratio ]
syms        = [ :wing, :htail, :vtail ]

other_horsies = [ atail_l_horsies[:]; atail_r_horsies[:] ]

# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = (
                                    wing    = data.wing.circulations, 
                                    htail   = data.htail.circulations, 
                                    vtail   = data.vtail.circulations, 
                                    atail_l = data.atail_l.circulations, 
                                    atail_r = data.atail_r.circulations
                                   ),
                    structures   = (
                                    wing  = Δx_wing, 
                                    htail = Δx_htail, 
                                    vtail = Δx_vtail
                                   ),
                    load_factor  = deg2rad(α))

# Set up initial guess and function
solve_aerostructural_residual!(R, x) =
    solve_coupled_residual!(R, x,
                            V, deg2rad(β), ρ, Ω,
                            syms, vlm_meshes, cam_meshes, fem_meshes,
                            other_horsies,
                            stiffy, weight, load_factor)

## Solve nonlinear system
reset_timer!()
@timeit "Solving Residuals" res_aerostruct =
    nlsolve(solve_aerostructural_residual!, x0,
            method         = :newton,
            show_trace     = true,
            # extended_trace = true,
            store_trace    = true,
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

new_horsies = new_horseshoes.(dxs, Ts, vlm_meshes, cam_meshes, fem_meshes)
all_horsies = [ reduce(vcat, vec.(new_horsies)); other_horsies ];

## Aerodynamic forces and center locations
U_opt      = freestream_to_cartesian(-V, α_opt, deg2rad(β))
new_acs    = new_horsies .|> horsies -> bound_leg_center.(horsies)
all_forces = nearfield_forces(Γ_opt, all_horsies, Γ_opt, all_horsies, U_opt, Ω, ρ)

new_Γs     = getindex.(Ref(Γ_opt), syms)
new_forces = getindex.(Ref(all_forces), syms)

## New beams and loads
new_fem_meshes = make_beam_mesh.(new_vlm_meshes, fem_weights)
fem_loads      = compute_loads.(new_acs, new_forces, new_fem_meshes);

## Compute stresses
δxs = Δs .|> dx -> eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = Δs .|> dx -> eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes
σs  = map((tubes, δx, δθ) -> reduce(hcat, von_mises_stress.(tubes, δx, δθ)), [ tubes_wing, tubes_htail, tubes_vtail ], δxs, δθs)

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
loads_plot   = fem_loads
σs_max       = map(σ -> maximum.(eachcol(σ)), σs)
σs_norm      = [ σ_max ./ maximum(σ_max) for σ_max in σs_max ]
σ_norms      = [ [ σ_norm; σ_norm[end] ] for σ_norm in σs_norm ]

## Panels
all_panels = ComponentArray(
                            wing    = wing_panels,
                            htail   = htail_panels,
                            vtail   = vtail_panels,
                            atail_l = atail_l_panels,
                            atail_r = atail_r_panels,
                           )

panel_plots = plot_panels(all_panels) 

# Aerodynamic centers and forces
ac_plots    = @. reduce(hcat, new_acs)
force_plots = @. reduce(hcat, new_forces)

# Cambers
cam_panels   = @. make_panels(cam_meshes)
cam_plot     = plot_panels(reduce(vcat, vec.(cam_panels)))
new_cam_plot = plot_panels(reduce(vcat, vec.(new_cam_panels)))

# Displacements
new_vlm_mesh_plot = @. reduce(hcat, new_vlm_meshes)
new_panel_plot    = plot_panels(reduce(vcat, vec.(make_panels.(new_vlm_meshes))))

## Planforms
wing_plan     = plot_wing(wing)
htail_plan    = plot_wing(htail)
vtail_plan    = plot_wing(vtail)
atail_l_plan  = plot_wing(atail_l)
atail_r_plan  = plot_wing(atail_r)

## New planforms
nwing_plan  = plot_wing(new_cam_meshes[1])
nhtail_plan = plot_wing(new_cam_meshes[2])
nvtail_plan = plot_wing(new_cam_meshes[3])

# Streamlines
seed    = [ 
            chop_coordinates(new_cam_meshes[1][1,:], 4)[1:2:end]; 
            # chop_coordinates(atail_l_cam_mesh[1,:], 4)[1:2:end] 
          ]
streams = plot_streams(fs, seed, all_horsies, Γ_opt, 5, 20);

## Plot
using Plots
using LaTeXStrings

# gr()
# /plotlyjs(size = (1280, 720))
pyplot(dpi = 150)
# pgfplotsx(size = (900, 600))

b = span(wing)
aircraft_plot =
    plot(xaxis = "x", yaxis = "y", zaxis = "z",
         camera = (45, 45),
        #  xlim = (-b/4, 3b/4),
     #     ylim = (-b/2, b/2),
         zlim = (-b/4, b/4),
         bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :topright,
         title = "Coupled Aerostructural Analysis"
        )

## Panels
[ plot!(pans, color = :lightgray, label = ifelse(i == 1, "Original Panels", :none)) for (i, pans) in enumerate(panel_plots) ]
[ plot!(pans, color = :lightblue, label = ifelse(i == 1, "Deflected Panels", :none)) for (i, pans) in enumerate(new_cam_plot) ]

## Undeflected planforms
plot!(wing_plan, color = :gray, label = "Original Wing")
plot!(htail_plan, color = :gray, label = "Horizontal Tail")
plot!(vtail_plan, color = :gray, label = "Vertical Tail")
plot!(atail_l_plan, color = :gray, label = "Left Taileron")
plot!(atail_r_plan, color = :gray, label = "Right Taileron")

# Deflected planforms
plot!(nwing_plan, color = :blue, label = "Deflected Wing")
plot!(nhtail_plan, color = :blue, label = "Deflected Horizontal Tail")
plot!(nvtail_plan, color = :blue, label = "Deflected Vertical Tail")

# Beams
thickness = 2.5
normer(rs) = [ rs; rs[end] ] / maximum(rs)
r_norms = @. normer([ r_wing, r_htail, r_vtail ]) * thickness
[ plot!(fem[1,:], fem[2,:], fem[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_ns) for (fem, r_ns) in zip(fem_plot, r_norms) ]

[ plot!(new_fem[1,:], new_fem[2,:], new_fem[3,:], m = (thickness, 0.8, :heat, Plots.stroke(0)), zcolor = σ_ns, cbar = true, label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_ns) for (new_fem, σ_ns, r_ns) in zip(new_fem_plot, σ_norms, r_norms) ]

# Streamlines
[ plot!(stream, color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(streams) ]

# Forces
# [ quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
#           quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.05,
#           label = ifelse(i == 1, "Panel Forces", :none), color = :orange) 
#           for (i, ac_plot, force_plot) in zip(1:length(ac_plots), ac_plots, force_plots) ]
# [ scatter!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], label = ifelse(i == 1, "Aerodynamic Centers", :none)) for ac_plot in ac_plots ]
# [ quiver!(fem_plots[1,:], fem_plots[2,:], fem_plots[3,:],
#         quiver=(loads_plots[1,:], loads_plots[2,:], loads_plots[3,:] ) .* 0.1,
#         label = ifelse(i == 1, "Beam Forces", :none)) for (fem_plots, loads_plots) in zip(fem_plot, loads_plot) ]

# savefig(aircraft_plot, "plots/AerostructWingTails.pdf")
plot!()