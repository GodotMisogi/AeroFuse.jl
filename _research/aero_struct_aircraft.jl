##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using Einsum

include("../src/Aerostructural/aerostruct.jl")

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [0.0, 0.0, 0.0],
            spans     = [5.0, 0.4],
            dihedrals = [0., 45.],
            sweep_LEs = [2.29, 60.]);

# Horizontal tail 
htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweep_LEs = [6.39])

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
                                         position = [4., 0, 0],
                                         angle    = deg2rad(-2.),
                                         axis     = [0., 1., 0.]
                                        )
vtail_panels, vtail_normals = panel_wing(vtail, 6, 4; 
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
aero_state.speed   = 30.
aero_state.alpha   = deg2rad(1.)
aero_state.beta    = deg2rad(0.)
aero_state.rho_ref = 0.98

# Build system with initial guess from aerodynamic-only analysis
aero_system = solve_case(aircraft, aero_state)
aero_surfs  = surfaces(aero_system)
print_coefficients(aero_surfs[1], aero_state);

horses     = horseshoes(aero_system)
Γ_0        = circulations(aero_system)

## Aerodynamic forces and center locations
vlm_forces = surface_forces(aero_surfs[1]) 
horsies    = horseshoes(aero_surfs[1])
vlm_acs    = bound_leg_center.(horsies)

## Mesh setup
vlm_mesh   = chord_coordinates(wing, span_num, chord_num)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 40 * 9.81
load_factor = 0.8;

## Structural variables

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3 # Density, kg/m³
ν     = 0.3   # Poisson ratio (UNUSED FOR NOW)
r     = 2e-2  # Outer radius, m
t     = 4e-3  # Thickness, m

# FEM mesh
fem_w    = 0.50
fem_mesh = make_beam_mesh(vlm_mesh, fem_w)
axes     = axis_transformation(fem_mesh, vlm_mesh)
Ls       = norm.(diff(fem_mesh)) # Beam lengths, m 
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, r, t)

# Stiffness matrix, loads and constraints
D         = build_big_stiffy(tubes, fem_mesh, vlm_mesh)
cons      = [length(fem_mesh) ÷ 2]
K         = build_stiffness_matrix(D, cons)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

## Solve system
dx        = solve_cantilever_beam(D, fem_loads, cons)
Δx        = [ zeros(6); dx[:] ]

## Aerostructural residual
#==========================================================================================#

solve_aerostructural_residual!(R, x) = solve_coupled_residual!(R, x, aero_system, aero_state, vlm_mesh, fem_mesh, K, weight, load_factor)

x0   = [ Γ_0; Δx; aero_state.alpha ]
res_aerostruct = nlsolve(solve_aerostructural_residual!, x0,
                         method     = :newton,
                         show_trace = true,
                        #  ftol       = 1e-7,
                        #  extended_trace = true,
                        #  autodiff   = :forward,
                        )

## Check numbers
lift     = total_force(values(aero_surfs))[3]
load_fac = lift / weight

println("Load factor: $load_fac")
println("Weight: $weight N")
println("Lift: $lift N")
println("Speed: $(aero_state.speed) m/s")
println("Angle of attack: $(rad2deg(aero_state.alpha))ᵒ")

## Get zero
x_opt = res_aerostruct.zero

# Get circulations
horsies = horseshoes(aero_system) 
Γs      = circulations(aero_system)

## Compute displacements
dx  = @views reshape(x_opt[length(Γ_0)+7:end-1], 6, length(fem_mesh))
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

## Generate DataFrame
df = DataFrame(permutedims([ fem_loads; dx ]), :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz, :dx, :dy, :dz, :θx, :θy, :θz])

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

# Beam circles
n_pts          = 20
circle3D(r, n) = [ (r*cos(θ), 0, r*sin(θ)) for θ in 0:2π/n:2π ];
circ           = circle3D(r * 5e-2, n_pts) 

draw_tube(p1, p2, circ) = [ [ circ_pt .+ p1, circ_pt .+ p2 ] for circ_pt in circ ]

beam_pts     = zip(tupvector(fem_mesh[1:end-1]), tupvector(fem_mesh[2:end]))
circ_pts     = [ draw_tube(pt[1], pt[2], circ) for pt in beam_pts ]

new_beam_pts = zip(tupvector(new_fem_mesh[1:end-1]), tupvector(new_fem_mesh[2:end]))
new_circ_pts = [ draw_tube(pt[1], pt[2], circ) for pt in new_beam_pts ]

# Beam loads
fem_plot   = reduce(hcat, fem_mesh)
loads_plot = fem_loads

# Aerodynamic centers and forces
panel_plot = plot_panels([ wing_panels[:]; htail_panels[:]; vtail_panels[:] ])
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

# Planforms and streams
wing_plan  = plot_wing(wing)
nwing_plan = plot_wing(new_vlm_mesh)
htail_plan = plot_wing(htail, 
                       position = [4., 0, 0],
                       angle    = deg2rad(-2.),
                       axis     = [0., 1., 0.]
                      )
vtail_plan = plot_wing(vtail, 
                       position = [4., 0, 0],
                       angle    = π/2, 
                       axis     = [1., 0., 0.]
                      )
fs         = Freestream(aero_state.speed, rad2deg(aero_state.alpha), rad2deg(aero_state.beta), aero_state.omega)
leading    = chop_coordinates(new_vlm_mesh[end,:], 1)
streams    = plot_streams(fs, leading, horsies, Γs, 5, 100)

## Plot
b = aero_state.span_ref
plot(
     camera = (-70, 20), 
     xlim = (-b/4, 3b/4),
    #  ylim = (-b/2, b/2), 
     zlim = (-b/8, b/4),
    )

# Panels
# [ plot!(pans, color = :gray, label = ifelse(i == 1, "Original", :none)) for (i, pans) in enumerate(panel_plot) ]
# [ plot!(pans, color = :blue, label = ifelse(i == 1, "Deflection", :none)) for (i, pans) in enumerate(new_panel_plot) ]

# Planforms
plot!(wing_plan, color = :black, label = "Original Wing")
plot!(nwing_plan, color = :blue, label = "Deflected Wing")
plot!(htail_plan, color = :brown, label = "Horizontal Tail")
plot!(vtail_plan, color = :brown, label = "Vertical Tail")

# Beams
[ plot!(reduce(vcat, pt), color = :gray, label = ifelse(i == 1, "Original Beam", :none)) for (i, pt) in enumerate(circ_pts) ]
[ plot!(reduce(vcat, pt), color = :cyan, label = ifelse(i == 1, "Deflected Beam", :none)) for (i, pt) in enumerate(new_circ_pts) ]

# Streamlines
[ plot!(stream,  color = :darkcyan, label = ifelse(i == 1, "Streamlines", :none)) for (i, stream) in enumerate(streams) ]

plot!()