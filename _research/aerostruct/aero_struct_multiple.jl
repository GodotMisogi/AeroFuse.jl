##

using AeroMDAO
using LinearAlgebra
using DataFrames
using NLsolve
using TimerOutputs
using SparseArrays
using ProfileView

# Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = Wing(foils     = Foil.(fill(naca4((2,4,1,2)), 2)),
            chords    = [1.0, 0.6],
            twists    = [0.0, 0.0],
            spans     = [5.0],
            dihedrals = [5.],
            LE_sweeps = [15.]);

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
vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.]);

## Meshing
wing_mesh  = WingMesh(wing, [6], 6);
htail_mesh = WingMesh(htail, [6], 3);
vtail_mesh = WingMesh(vtail, [4], 3);

## Aircraft assembly
aircraft = ComponentVector(
                           wing  = make_horseshoes(wing_mesh),
                           htail = make_horseshoes(htail_mesh),
                           vtail = make_horseshoes(vtail_mesh),
                          );

## Aerodynamic case

# Freestream conditions
fs      = Freestream(alpha = 3.0, 
                     beta  = 0.0, 
                     omega = [0.,0.,0.])

# Reference values
refs    = References(
                     speed    = 25.0,
                     density  = 0.98, 
                     area     = projected_area(wing),   
                     span     = span(wing), 
                     chord    = mean_aerodynamic_chord(wing), 
                     location = mean_aerodynamic_center(wing)
                    )

## Solve aerodynamic case for initial vector
@time system = solve_case(aircraft, fs, refs;
                        #   print_components = true,
                         )

## Data collection
@time Fs = surface_forces(system);

## Wing FEM setup
vlm_acs_wing    = bound_leg_center.(system.vortices.wing)
vlm_forces_wing = Fs.wing

wing_beam_ratio = 0.40
wing_fem_mesh   = make_beam_mesh(wing_mesh.chord_mesh, wing_beam_ratio)

aluminum = Material(       # Aluminum properties
                    85e9,  # Elastic modulus, N/m²
                    25e9,  # Shear modulus, N/m²,
                    350e6, # Yield stress with factor of safety 2.5, N/m²,
                    1.6e3, # Density, kg/m³
                    )

Ls_wing = norm.(diff(wing_fem_mesh))                 # Beam lengths, m
rs_wing = LinRange(2e-2, 1e-2, length(Ls_wing) ÷ 2)  # Outer radius, m
ts_wing = LinRange(1e-2, 6e-3, length(Ls_wing) ÷ 2)  # Thickness, m
r_wing  = [ reverse(rs_wing); rs_wing ]
t_wing  = [ reverse(ts_wing); ts_wing ]

tubes_wing     = Tube.(Ref(aluminum), Ls_wing, r_wing, t_wing)
wing_beam      = Beam(tubes_wing)


as_wing = AerostructWing(wing_mesh, wing_beam)

Ks_wing        = build_big_stiffy(tubes_wing, wing_fem_mesh, wing_mesh.chord_mesh)
cons_wing      = [length(wing_fem_mesh) ÷ 2]
stiffy_wing    = build_stiffness_matrix(Ks_wing, cons_wing)
fem_loads_wing = compute_loads(vlm_acs_wing, vlm_forces_wing, wing_fem_mesh)

dx_wing = solve_cantilever_beam(Ks_wing, fem_loads_wing, cons_wing)
Δx_wing = [ zeros(6); dx_wing[:] ]

## Horizontal tail FEM setup
vlm_acs_htail    = bound_leg_center.(system.vortices.wing)
vlm_forces_htail = Fs.htail

htail_beam_ratio = 0.35
htail_fem_mesh   = make_beam_mesh(htail_mesh.chord_mesh, htail_beam_ratio)

# Beam properties
Ls_htail = norm.(diff(htail_fem_mesh))                              # Beam lengths, m
rs_htail = LinRange(8e-3, 4e-3, length(Ls_htail) ÷ 2)  # Outer radius, m
ts_htail = LinRange(8e-4, 6e-4, length(Ls_htail) ÷ 2)  # Thickness, m
r_htail  = [ reverse(rs_htail); rs_htail ]
t_htail  = [ reverse(ts_htail); ts_htail ]

tubes_htail     = Tube.(Ref(aluminum), Ls_htail, r_htail, t_htail)
htail_beam      = Beam(tubes_htail)

Ks_htail        = build_big_stiffy(tubes_htail, htail_fem_mesh, htail_mesh.chord_mesh)
cons_htail      = [length(htail_fem_mesh) ÷ 2]
stiffy_htail    = build_stiffness_matrix(Ks_htail, cons_htail)
fem_loads_htail = compute_loads(vlm_acs_htail, vlm_forces_htail, htail_fem_mesh)

dx_htail = solve_cantilever_beam(Ks_htail, fem_loads_htail, cons_htail)
Δx_htail = [ zeros(6); dx_htail[:] ]

as_htail = AerostructWing(htail_mesh, htail_beam);

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

surfs       = [as_wing, as_htail] 

chord_meshes  = [ wing_mesh.chord_mesh, htail_mesh.chord_mesh ]
camber_meshes  = [ wing_mesh.camber_mesh, htail_mesh.camber_mesh ]
fem_meshes  = [ wing_fem_mesh, htail_fem_mesh ]
fem_weights = [ wing_beam_ratio, htail_beam_ratio ]
syms        = [ :wing, :htail ]

# Initial guess as ComponentArray for the different equations
x0 = ComponentArray(aerodynamics = (
                                    wing  = system.circulations.wing, 
                                    htail = system.circulations.htail, 
                                    vtail = system.circulations.vtail
                                   ),
                    structures   = (
                                    wing  = Δx_wing, 
                                    htail = Δx_htail
                                   ),
                    load_factor  = deg2rad(α))

# Set up initial guess and function
function aerostructural_problem(V, β, ρ, Ω, syms, chord_meshes, camber_meshes, fem_meshes, horsies, stiffy, weight, load_factor)
    f!(R, x) =
        solve_coupled_residual!(R, x,
                                V, β, ρ, Ω,         # Aerodynamic state
                                syms, chord_meshes, camber_meshes, fem_meshes,
                                horsies, stiffy, weight, load_factor)
end

@code_warntype aerostructural_problem(refs.speed, fs.beta, refs.density, fs.omega, syms, chord_meshes, camber_meshes, fem_meshes, aircraft.vtail, stiffy, weight, load_factor)

## Closure
solve_aerostruct! = aerostructural_problem(refs.speed, fs.beta, refs.density, fs.omega, syms, chord_meshes, camber_meshes, fem_meshes, aircraft.vtail, stiffy, weight, load_factor)

## Solve nonlinear system
using ForwardDiff, ReverseDiff
# using Zygote

function newton_raphson(f!, x0; max_iters = 50, tol = 1e-9)
    x = copy(x0)
    dx = similar(x)
    R = similar(x)
    ∂R∂x = Matrix{eltype(x)}(undef, length(R), length(x))
    ε = 1e5
    i = 0
    for i = 1:max_iters
        @timeit "Evaluating Jacobian" ForwardDiff.jacobian!(∂R∂x, f!, R, x)
        @timeit "Inverting System" ldiv!(dx, factorize(∂R∂x), -R)
        if ε <= tol return x end # Needs NAN checks and everything like NLsolve
        ε    = LinearAlgebra.norm(dx)
        # @show (i, ε)
        x  .+= dx
        i   += 1
    end
    return x
end

##
R = similar(x0)
@code_warntype solve_aerostruct!(R, x0)

##
reset_timer!()

@timeit "Solving Residuals" res_aerostruct =
    # @time newton_raphson(solve_aerostruct!, x0)
    nlsolve(solve_aerostruct!, x0,
            method         = :newton,
            autodiff       = :forward,
            show_trace     = true,
            # extended_trace = true,
            )
# 
print_timer()

## Get zero
x_opt = res_aerostruct.zero
Γ_opt = x_opt.aerodynamics
δ_opt = x_opt.structures
α_opt = x_opt.load_factor

## Compute displacements
Δs    = map((key, n) -> reshape(δ_opt[key][7:end], 6, n), valkeys(δ_opt), length.(fem_meshes))
dxs   = mesh_translation.(Δs)
Ts    = mesh_rotation.(Δs)

## New VLM variables
new.chord_meshes = transfer_displacements.(dxs, Ts, chord_meshes, fem_meshes)
new_panels     = make_panels.(new.chord_meshes)

new_camber_meshes = transfer_displacements.(dxs, Ts, camber_meshes, fem_meshes)
new_cam_panels = make_panels.(new_camber_meshes)

new_horsies = new_horseshoes.(dxs, Ts, chord_meshes, camber_meshes, fem_meshes)
all_horsies = [ reduce(vcat, vec.(new_horsies)); other_horsies ];

## Aerodynamic forces and center locations
U_opt      = freestream_to_cartesian(-V, α_opt, deg2rad(β))
new_acs    = new_horsies .|> horsies -> bound_leg_center.(horsies)
all_forces = surface_forces(Γ_opt, all_horsies, Γ_opt, all_horsies, U_opt, Ω, ρ)

new_Γs     = getindex.(Ref(Γ_opt), syms)
new_forces = surface_forces.(new_Γs, new_horsies, Ref(Γs), Ref(all_horsies), Ref(U_opt), Ref(Ω), Ref(ρ));

## New beams and loads
new_fem_meshes = make_beam_mesh.(new.chord_meshes, fem_weights)
fem_loads      = compute_loads.(new_acs, new_forces, new_fem_meshes);

## Compute stresses
δxs = Δs .|> dx -> eachcol(diff(dx[1:3,:], dims = 2)) # Need to transform these to principal axes
δθs = Δs .|> dx -> eachcol(diff(dx[4:6,:], dims = 2)) # Need to transform these to principal axes
σs  = map((tubes, δx, δθ) -> reduce(hcat, von_mises_stress.(tubes, δx, δθ)), [ tubes_wing, tubes_htail ], δxs, δθs)

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
wing_panel_plot  = plot_panels(wing_panels)
htail_panel_plot = plot_panels(htail_panels)
vtail_panel_plot = plot_panels(vtail_panels)

# Aerodynamic centers and forces
ac_plot    = @. reduce(hcat, new_acs)
force_plot = @. reduce(hcat, new_forces)

# Cambers
cam_panels   = @. make_panels(camber_meshes)
cam_plot     = plot_panels(reduce(vcat, vec.(cam_panels)))
new_cam_plot = plot_panels(reduce(vcat, vec.(new_cam_panels)))

# Displacements
new.chord_mesh_plot = @. reduce(hcat, new.chord_meshes)
new_panel_plot    = plot_panels(reduce(vcat, vec.(make_panels.(new.chord_meshes))))

# Planforms
wing_plan   = plot_wing(wing)
htail_plan  = plot_wing(htail)
vtail_plan  = plot_wing(vtail)

# New planforms
nwing_plan  = plot_wing(new_camber_meshes[1])
nhtail_plan = plot_wing(new_camber_meshes[2])

# Streamlines
seed    = chop_coordinates(new_camber_meshes[1][end,:], 4)
streams = streamlines(fs, seed, all_horsies, Γ_opt, 5, 100);

## Plot
using Plots
using LaTeXStrings

# gr()
# plotlyjs(dpi = 300, size = (1280, 720))
pyplot(dpi = 300)
# pgfplotsx(size = (900, 600))

aircraft_plot =
    plot(xaxis = "x", yaxis = "y", zaxis = "z",
         camera = (-75, 20),
         xlim = (-b/4, 3b/4),
     #     ylim = (-b/2, b/2),
         zlim = (-b/8, b/4),
        #  bg_inside = RGBA(0.96, 0.96, 0.96, 1.0),
         legend = :topright,
         title = "Coupled Aerostructural Analysis"
        )

# Panels
[ plot!(pans, color = :lightgray, label = ifelse(i == 1, "Original Panels", :none), linestyle = :solid) for (i, pans) in enumerate(cam_plot) ]
[ plot!(pans, color = :lightblue, label = ifelse(i == 1, "Deflected Panels", :none), linestyle = :solid) for (i, pans) in enumerate(new_cam_plot) ]
[ plot!(pans, color = :brown, label = :none, linestyle = :solid) for (i, pans) in enumerate(vtail_panel_plot) ]

# Planforms
plot!(wing_plan, color = :gray, label = "Original Wing", linestyle = :solid)
plot!(nwing_plan, color = :blue, label = "Deflected Wing")
plot!(htail_plan, color = :gray, label = "Horizontal Tail")
plot!(nhtail_plan, color = :blue, label = "Deflected Horizontal Tail")
plot!(vtail_plan, color = :brown, label = "Vertical Tail")

# Beams
thickness = 2.5
normer(rs) = [ rs; rs[end] ] / maximum(rs)
r_norms = @. normer([ r_wing, r_htail ]) * thickness
[ plot!(fem[1,:], fem[2,:], fem[3,:], color = :black, label = "Original Beam", linestyle = :solid, linewidth = r_ns) for (fem, r_ns) in zip(fem_plot, r_norms) ]

[ plot!(new_fem[1,:], new_fem[2,:], new_fem[3,:], m = (thickness, 0.8, :heat, Plots.stroke(0)), zcolor = σ_ns, cbar = true, label = "Deflected Beam Stresses", linestyle = :solid, linewidth = r_ns) for (new_fem, σ_ns, r_ns) in zip(new_fem_plot, σ_norms, r_norms) ]

# Streamlines
[ plot!(stream, color = RGBA(0.5, 0.8, 0.5, 1.0), label = ifelse(i == 1, "Streamlines", :none), linestyle = :solid) for (i, stream) in enumerate(streams) ]

# Forces
# quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
#         quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1,
#         label = "Panel Forces", color = :orange)
# scatter!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], label = "Aerodynamic Centers")
# quiver!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:],
#         quiver=(loads_plot[1,:], loads_plot[2,:], loads_plot[3,:] ) .* 0.1,
#         label = "Beam Forces")

# savefig(aircraft_plot, "plots/AerostructWingTail.html")
plot!()