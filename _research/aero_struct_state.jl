##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using Einsum
using BenchmarkTools
using TimerOutputs

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
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
             sweep_LEs = [6.39])
htail_position = [4., 0, 0.]
htail_angle    = deg2rad(-1.)

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)), 
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97]);

## Meshing and assembly
span_num, chord_num         = [8, 3], 6
wing_panels,  wing_normals  = panel_wing(wing, span_num, chord_num;)
htail_panels, htail_normals = panel_wing(htail, 6, 6;
                                         position = htail_position,
                                         angle    = htail_angle,
                                         axis     = [0., 1., 0.]
                                        )
vtail_panels, vtail_normals = panel_wing(vtail, 6, 6; 
                                         position = [4., 0, 0],
                                         angle    = π/2, 
                                         axis     = [1., 0., 0.]
                                        )

# Aircraft assembly
aircraft = Dict(
                "Wing"            => (wing_panels,  wing_normals),
                "Horizontal Tail" => (htail_panels, htail_normals),
                "Vertical Tail"   => (vtail_panels, vtail_normals),
               );
wing_mac = mean_aerodynamic_center(wing);

# Set up aerodynamic state
aero_state = VLMState(0., 0., 0., [0.0, 0.0, 0.0], 
                      rho_ref   = 1.225,
                      r_ref     = [ wing_mac[1], 0., 0. ],
                      area_ref  = projected_area(wing), 
                      chord_ref = mean_aerodynamic_chord(wing), 
                      span_ref  = span(wing));

# Test case - Fixed speed
aero_state.speed   = 25.
aero_state.alpha   = deg2rad(3.)
aero_state.beta    = deg2rad(0.)
aero_state.rho_ref = 0.98

# Build system with initial guess from aerodynamic-only analysis
aero_system = solve_case(aircraft, aero_state)
aero_surfs  = surfaces(aero_system)
print_coefficients(aero_surfs[1], aero_state);

# Get initial aerodynamic vector for Newton method
Γ_0        = circulations(aero_system)

## Aerodynamic forces and center locations
horsies    = horseshoes(aero_surfs[1])
vlm_acs    = bound_leg_center.(horsies)
vlm_forces = surface_forces(aero_surfs[1]) 

## Mesh setup
vlm_mesh   = chord_coordinates(wing, span_num, chord_num)

# FEM mesh
fem_w    = 0.40
fem_mesh = make_beam_mesh(vlm_mesh, fem_w)
axes     = axis_transformation(fem_mesh, vlm_mesh)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 60 * 9.81
load_factor = 1.3;

## Structural variables

# Beam properties
Ls    = norm.(diff(fem_mesh)) # Beam lengths, m 
E     = 85e9                  # Elastic modulus, N/m²
G     = 25e9                  # Shear modulus, N/m²
σ_max = 350e6                 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3                 # Density, kg/m³
ν     = 0.3                   # Poisson ratio (UNUSED FOR NOW)
rs    = range(2e-2, stop = 8e-3, length = length(Ls) ÷ 2)  # Outer radius, m
ts    = range(8e-3, stop = 2e-3, length = length(Ls) ÷ 2)  # Thickness, m
r     = [ reverse(rs); rs ]
t     = [ reverse(ts); ts ]

aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, r, t)

# Stiffness matrix, loads and constraints
D         = build_big_stiffy(tubes, fem_mesh, vlm_mesh)
cons      = [length(fem_mesh) ÷ 2]
stiffy    = build_stiffness_matrix(D, cons)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

## Solve system and get initial structural vector for Newton method
dx        = solve_cantilever_beam(D, fem_loads, cons)
Δx        = [ zeros(6); dx[:] ]

## Aerostructural residual
#==========================================================================================#

# Set up initial guess and function
solve_aerostructural_residual!(R, x) = solve_coupled_residual!(R, x, aero_system, aero_state, vlm_mesh, fem_mesh, stiffy, weight, load_factor)

# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = Γ_0, structures = Δx, load_factor = deg2rad(α))

## Solve system
reset_timer!()
@time res_aerostruct = nlsolve(solve_aerostructural_residual!, x0,
                         method     = :newton,
                         show_trace = true,
                        #  ftol       = 1e-7,
                        #  extended_trace = true,
                        #  autodiff   = :forward,
                        );
print_timer()

## Check numbers
lift     = total_force(system, aero_state)[3]
load_fac = lift * cos(aero_state.alpha) / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(aero_state.speed) m/s")
println("Angle of attack: $(rad2deg(aero_state.alpha))ᵒ")

## Get zero
x_opt = res_aerostruct.zero
Γ_opt = x_opt.aerodynamics
δ_opt = x_opt.structures[7:end]
α_opt = x_opt.load_factor

# Get circulations
horsies = horseshoes(aero_system) 
Γs      = circulations(aero_system)

## Compute displacements
dx  = @views reshape(δ_opt, 6, length(fem_mesh))
dxs = @views SVector.(dx[1,:], dx[2,:], dx[3,:])
Ts  = rotation_matrix(dx[4:6,:])

# Perturb VLM mesh and normals
new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
new_panels   = make_panels(new_vlm_mesh)

# New beams
new_fem_mesh = make_beam_mesh(new_vlm_mesh, fem_w)

## Aerodynamic forces and center locations
horsies    = horseshoes(aero_surfs[1])
vlm_acs    = bound_leg_center.(horsies)
vlm_forces = surface_forces(aero_surfs[1])
fem_loads  = compute_loads(vlm_acs, vlm_forces, new_fem_mesh)

## Compute stresses
δxs = eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes 
σs  = reduce(hcat, von_mises_stress.(tubes, δxs, δθs))

## Generate DataFrame
df = DataFrame(permutedims([ fem_loads; dx ]), :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz, :dx, :dy, :dz, :θx, :θy, :θz])

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

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

# Displacements
new_vlm_mesh_plot = reduce(hcat, new_vlm_mesh)
new_panel_plot = plot_panels(make_panels(new_vlm_mesh)[:])

xs_plot = reduce(hcat, (fem_mesh[1:end-1] + fem_mesh[2:end]) / 2)
axes    = axis_transformation(fem_mesh, vlm_mesh)
cs_plot = axes[:,1,:]
ss_plot = axes[:,2,:]
ns_plot = axes[:,3,:]

# Planforms
wing_plan  = plot_wing(wing)
nwing_plan = plot_wing(new_vlm_mesh)
htail_plan = plot_wing(htail, 
                       position = htail_position,
                       angle    = htail_angle,
                       axis     = [0., 1., 0.]
                      )
vtail_plan = plot_wing(vtail, 
                       position = [4., 0, 0],
                       angle    = π/2, 
                       axis     = [1., 0., 0.]
                      )

# Streamlines
fs      = Freestream(aero_state.speed, rad2deg(aero_state.alpha), rad2deg(aero_state.beta), aero_state.omega)
seed    = chop_coordinates(new_vlm_mesh[end,:], 2)
streams = plot_streams(fs, seed, horsies, Γs, 2.5, 100);

## Plot
b = aero_state.span_ref
using LaTeXStrings
# pgfplotsx(size = (900, 600))
aircraft_plot = 
    plot(xaxis = L"$x$", yaxis = L"$y$", zaxis = L"$z$",
         camera = (-75, 30), 
         xlim = (0, b/2),
     #     ylim = (-b/2, b/2),
         zlim = (-b/8, b/8),
         bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :bottomright,
         title = "Coupled Aerostructural Analysis"
        )

# Panels
# [ plot!(pans, color = :gray,  label = ifelse(i == 1, "Original Wing Panels", :none),  linestyle = :solid) for (i, pans) in enumerate(wing_panel_plot)  ]
# [ plot!(pans, color = RGBA(0.5, 0.5, 0.8, 0.7),  label = ifelse(i == 1, "Deflected Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_panel_plot)   ]
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
plot!(new_fem_plot[1,:], new_fem_plot[2,:], new_fem_plot[3,:], color = RGBA.(σ_norms, 0.5, 0.6, 1.0), label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_norm)

# Streamlines
[ plot!(stream,  color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(streams) ]

# Forces
# quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
#         quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1,
#         label = "Panel Forces", color = :orange)

# savefig(aircraft_plot, "plots/AerostructAircraft.pdf")
plot!()