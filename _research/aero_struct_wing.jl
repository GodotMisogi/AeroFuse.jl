##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve
using Einsum
using TimerOutputs

include("../src/Aerostructural/aerostruct.jl")

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = WingSection(root_foil  = naca4((0,0,1,2)),
                   span       = 2.6,
                   dihedral   = 1.0,
                   sweep_LE   = 20.0,
                   taper      = 0.5,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
print_info(wing)

# Mesh details
span_num  = [12]
chord_num = 6

## Aerodynamic case
ρ       = 0.98
ref     = [0.25 * wing_mac[1], 0., 0.]
V, α, β = 25.0, 5.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

# Solve aerodynamic case for initial vector
@time nf_coeffs, ff_coeffs, CFs, CMs, panels, normies, horsies, Γs = 
    solve_case(wing, fs;
               rho_ref   = ρ,
               r_ref     = ref,
               span_num  = span_num[1] * 2,
               chord_num = chord_num,
               spacing   = ["sine"]
              );

print_coefficients(nf_coeffs, ff_coeffs)

## Aerodynamic forces and center locations
vlm_acs    = bound_leg_center.(horsies)
vlm_forces = force.(CFs, dynamic_pressure(ρ, V), projected_area(wing))

## Mesh setup
vlm_mesh   = chord_coordinates(wing, span_num, chord_num, spacings = ["sine"])

# FEM mesh
fem_w    = 0.40
fem_mesh = make_beam_mesh(vlm_mesh, fem_w)
axes     = axis_transformation(fem_mesh, vlm_mesh)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 10 * 9.81
load_factor = 1.0;

## Structural variables

## Beam properties
Ls    = norm.(diff(fem_mesh)) # Beam lengths, m 
E     = 85e9                  # Elastic modulus, N/m²
G     = 25e9                  # Shear modulus, N/m²
σ_max = 350e6                 # Yield stress with factor of safety 2.5, N/m²
rho   = 1.6e3                 # Density, kg/m³
ν     = 0.3                   # Poisson ratio (UNUSED FOR NOW)
rs    = range(6e-3, stop = 4e-3, length = length(Ls) ÷ 2)  # Outer radius, m
ts    = range(1e-3, stop = 8e-4, length = length(Ls) ÷ 2)  # Thickness, m
r     = [ reverse(rs); rs ]
t     = [ reverse(ts); ts ]

aluminum = Material(E, G, σ_max, rho)
tubes    = Tube.(Ref(aluminum), Ls, r, t)

## Stiffness matrix, loads and constraints
D         = build_big_stiffy(tubes, fem_mesh, vlm_mesh)
cons      = [length(fem_mesh) ÷ 2]
stiffy    = build_stiffness_matrix(D, cons)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

## Solve system
dx = solve_cantilever_beam(D, fem_loads, cons)
Δx = [ zeros(6); dx[:] ]

## Aerostructural residual
#==========================================================================================#

solve_aerostructural_residual!(R, x) = 
    solve_coupled_residual!(R, x, 
                            V, deg2rad(β), ρ, Ω, # Aerodynamic state
                            vlm_mesh, normies,   # Aerodynamic variables
                            fem_mesh, stiffy,    # Structural variables
                            weight, load_factor) # Load factor variables

x0   = [ Γs[:]; Δx; deg2rad(α) ]
@time res_aerostruct = nlsolve(solve_aerostructural_residual!, x0,
                         method     = :newton,
                         show_trace = true,
                        #  extended_trace = true,
                         autodiff   = :forward,
                        )

## Get zero
x_opt = res_aerostruct.zero
Γ_opt = @views reshape(x_opt[1:length(normies)], size(normies))
δ_opt = @view x_opt[length(normies)+7:end-1]
α_opt = x_opt[end];

## Compute displacements
dx  = @views reshape(δ_opt, 6, length(fem_mesh))
dxs = @views SVector.(dx[1,:], dx[2,:], dx[3,:])
Ts  = rotation_matrix(dx[4:6,:])

## Perturb VLM mesh and normals
new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
new_panels   = make_panels(new_vlm_mesh)

## Aerodynamic forces and center locations
U_opt       = freestream_to_cartesian(-V, α_opt, deg2rad(β))
new_horsies = horseshoe_line.(new_panels)
vlm_acs     = bound_leg_center.(new_horsies)
vlm_forces  = nearfield_forces(Γ_opt, new_horsies, Γ_opt, new_horsies, U_opt, Ω, ρ)

## New beams and loads
new_fem_mesh = make_beam_mesh(new_vlm_mesh, fem_w)
fem_loads    = compute_loads(vlm_acs, vlm_forces, new_fem_mesh)

## Compute stresses
δxs = eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes 
σs  = reduce(hcat, von_mises_stress.(tubes, δxs, δθs))

## Check numbers
lift     = sum(vlm_forces)[3]
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

using Plots
pyplot(dpi = 300)

# Beam loads and stresses
fem_plot     = reduce(hcat, chop_coordinates(fem_mesh, 1))
new_fem_plot = reduce(hcat, chop_coordinates(new_fem_mesh, 1))
loads_plot   = fem_loads
σs_max       = maximum.(eachcol(σs))
σs_norm      = σs_max ./ maximum(σs_max)
σ_norms      = [ σs_norm; σs_norm[end] ]

# Aerodynamic centers and forces
panel_plot = plot_panels(panels[:])
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

# Streamlines
seed    = chop_coordinates(new_vlm_mesh[end,:], 2)
streams = plot_streams(fs, seed, new_horsies, Γs, 2.5, 100);

## Plot
using LaTeXStrings
pgfplotsx(size = (900, 600))
aircraft_plot = 
    plot(xaxis = L"$x$", yaxis = L"$y$", zaxis = L"$z$",
         camera = (-80, 30), 
         xlim = (-b/4, 3b/4),
     #     ylim = (-b/2, b/2),
         zlim = (-b/8, b/4),
         bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :bottomright,
         title = "Coupled Aerostructural Analysis"
        )


# Panels
# [ plot!(pans, color = :gray,  label = ifelse(i == 1, "Original Wing Panels", :none),  linestyle = :solid) for (i, pans) in enumerate(panel_plot)  ]
[ plot!(pans, color = RGBA(0.5, 0.5, 0.8, 0.7),  label = ifelse(i == 1, "Deflected Wing Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_panel_plot)   ]

# Planform
plot!(wing_plan,  color = :gray, label = "Original Wing", linestyle = :solid)
plot!(nwing_plan, color = :blue, label = "Deflected Wing")

# Beams
thickness = 2.5
r_norm = [ r; r[end]] / maximum(r) * thickness
plot!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_norm)
plot!(new_fem_plot[1,:], new_fem_plot[2,:], new_fem_plot[3,:], color = RGBA.(σ_norms, 0.5, 0.6, 1.0), label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_norm)

# Streamlines
[ plot!(stream,  color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(streams) ]

# Forces
quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
        quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1,
        label = "Panel Forces", color = :orange)
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

savefig(aircraft_plot, "plots/AerostructWing.pdf")
# plot!()