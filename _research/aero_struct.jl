##
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using NLsolve

## Helper functions
#==========================================================================================#

# Sum adjacent values
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

## Aerodynamic analysis
#==========================================================================================#

# lifting_line_forces(surf) = sum(surface_forces(surf), dims = 1)[:]
aerodynamic_forces(surfs) = reduce(vcat, (vec ∘ surface_forces).(surfs))

points(horses) = @. r1(bound_leg_center(horses), horses)[:], r2(bound_leg_center(horses), horses)[:]

total_force(surfs) = sum(sum ∘ surface_forces, surfs) 

function compute_loads(forces, r1s, r2s)
    # Aerodynamic forces
    half_forces = @. forces / 2
    M1s         = @. r1s × half_forces
    M2s         = @. r2s × half_forces

    # Interspersing
    forces  = adjacent_joiner(half_forces, half_forces)
    moments = adjacent_joiner(M1s, M2s)

    # Boundary condition, setting F = 0 at center of wing
    n = ceil(Int, length(horses) / 2) + 1
    forces[n-1:n+1] .-= forces[n-1:n+1];

    # Assembly
    Fx = getindex.(forces, 1)
    Fy = getindex.(forces, 2)
    Fz = getindex.(forces, 3)
    Mx = getindex.(moments, 1)
    My = getindex.(moments, 2)
    Mz = getindex.(moments, 3)

    py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
    px = [ Fx[:]; Mx[:] ]

    [ py; pz; px ]
end

## Structural analysis
#==========================================================================================#

solve_beam_residual!(R, K, δ, F) = R .= K * δ - F

function compute_displacements(Δs, r1s, r2s)
    # Assembly
    n  = length(r1s) + 1    # Number of finite element points
    n1 = 2n                 # dim(Fy + My)
    n2 = n1 + 2n            # dim(Fz + Mz)
    n3 = n2 +  n            # dim(Fx)
    n4 = n3 +  n            # dim(Mx)

    dy = @view Δs[1:2:n1]
    θy = @view Δs[2:2:n1]
    dz = @view Δs[n1+1:2:n2]
    θz = @view Δs[n1+2:2:n2]
    dx = @view Δs[n2+1:n3]
    θx = @view Δs[n3+1:n4]

    # Coordinates
    rs = @. SVector(dx, dy, dz)
    θs = @. SVector(θx, θy, θz)
    ds = rs .+ θs .× adjacent_joiner(r1s, r2s)

    # A bijective mapping between the wing's geometric and coordinate representations needs to be defined, if it exists.
    # Suspecting it doesn't due to different perturbations in the section spans at the leading and trailing edges.

    ds, θs
end

## Load-displacement transfer mechanisms
#==========================================================================================#

function transfer_displacements!(system :: VLMSystem, xyzs, normals, ds, θs)
    # Displace and rotate about quarter-chord point, where the beam is located. Needs many corrections...
    Δxyzs      = permutedims([ ds ds ])
    new_xyzs   = xyzs + Δxyzs
    new_pans   = make_panels(new_xyzs)
    new_horses = horseshoe_line.(new_pans) 
    # new_norms  = panel_normal.(new_pans)
    new_norms  = normals[:] 
    # + (θs[1:end-1] + θs[2:end]) / 2 # WRONG

    # Set new system variables
    system.horseshoes = new_horses[:]
    system.normals    = new_norms[:]

    new_xyzs
end

## Weights and fuel loads
#==========================================================================================#

load_factor_residual(L, W, n) = L - n * W
load_factor_residual!(R, L, W, n) = R .= load_factor_residual(L, W, n)

## Coupled residual system
#==========================================================================================#

function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_surfs :: Vector{<: VLMSurface}, aero_state :: VLMState, xyzs, stiffness_matrix, weight, load_factor)
    n = (prod ∘ size)(horseshoes(aero_system))   # Get panel size

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute displacements
    r1s, r2s = (points ∘ horseshoes)(aero_system)
    ds, θs   = compute_displacements(δ, r1s, r2s)

    # The displacements needs to be transformed into the coordinate frame before transferring...
    
    # Panel local coordinate system
    horsies = horseshoes(aero_system)
    bounds  = bound_leg_vector.(horsies)
    normies = Ref(aero_state.velocity) .× bounds

    # Compute displacements
    transfer_displacements!(aero_system, xyzs, normies, ds, θs)

    # Aerodynamic residuals
    aero_state.alpha = x[end]
    solve_aerodynamic_residual!(R_A, Γ, aero_system, aero_surfs, aero_state)
    
    # Compute forces
    forces = aerodynamic_forces(aero_surfs)

    # The forces need to be transformed into the principal axes of the finite-elements...
    # Also need to consider boundary conditions...

    # Direction cosine transformations
    streams = @. bounds × normies
    glob    = reduce(hcat, fill([[1 0 0], [0 1 0], [0 0 1]], 3))
    dircos  = [ dot.(reduce(hcat, fill([c', n', s'], 3)), glob)  for (c, n, s) in zip(bounds, normies, streams) ]
    F_S     = dircos .* forces
    F       = compute_loads(F_S, r1s, r2s)

    # Compute lift for load factor residual
    L = total_force(aero_surfs)[3]

    # Structural residuals
    solve_beam_residual!(R_S, stiffness_matrix, δ, F)

    # Weight residual
    load_factor_residual!(R_W, L, weight, load_factor)

    R
end

## Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = WingSection(root_foil  = naca4((0,0,1,2)),
                   span       = 1.3,
                   dihedral   = 5.0,
                   sweep_LE   = 0.0,
                   taper      = 1.0,
                   root_chord = 0.314,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
wing_name   = "Wing"
print_info(wing)

# Mesh
span_num        = 5
chord_num       = 1
xyzs            = chord_coordinates(wing, [span_num], chord_num)
panels, normies = panel_wing(wing, span_num, chord_num);
aircraft        = Dict(wing_name => (panels, normies));

# Set up aerodynamic state
aero_state = VLMState(0., 0., 0., [0.0, 0.0, 0.0], 
                      rho_ref   = 1.225,
                      r_ref     = [ wing_mac[1], 0., 0. ],
                      area_ref  = projected_area(wing), 
                      chord_ref = mean_aerodynamic_chord(wing), 
                      span_ref  = span(wing));

# Test case - Fixed speed
aero_state.speed       = 20.
aero_state.alpha   = deg2rad(1.)
aero_state.rho_ref = 0.98

# Build system with initial guess from aerodynamic-only analysis
aero_system, aero_surfs = solve_case(aircraft, aero_state)
print_coefficients(aero_surfs[1], aero_state);

horses  = horseshoes(aero_system)
Γ_0     = circulations(aero_system)
normies = normals(aero_system)

## Weight variables (FOR FUTURE USE)

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 15 * 9.81
load_factor = 1.0;

## Structural variables

# Material - Aluminum
E     = 70e9  # Elastic modulus, N/m²
G     = 30e9  # Shear modulus, N/m²
σ_max = 200e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 3e3   # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-1  # Outer radius, m
t     = 1e-3  # Thickness, m

Ls    = @. (norm ∘ bound_leg_vector)(horses) # Length, m

aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

## Stiffness matrix setup
K        = tube_stiffness_matrix(aluminum, tubes)
forces   = aerodynamic_forces(aero_surfs)
r1s, r2s = (points ∘ horseshoes)(aero_system)

# Direction cosine matrix transformation
streams  = @. bounds × normies
glob     = reduce(hcat, fill([[1 0 0], [0 1 0], [0 0 1]], 3))
dircos   = [ dot.(reduce(hcat, fill([c', n', s'], 3)), glob)  for (c, n, s) in zip(bounds, normies, streams) ]
F_S      = dircos .* forces
F        = compute_loads(F_S, r1s, r2s)

# Generate initial guess from structural-only analysis
Δx       = K \ F

## Aerostructural residual
#==========================================================================================#

solve_aerostructural_residual!(R, x) = solve_coupled_residual!(R, x, aero_system, (collect ∘ values)(aero_surfs), aero_state, xyzs, K, weight, load_factor)

x0   = [ Γ_0; Δx; aero_state.alpha ]
res_aerostruct = nlsolve(solve_aerostructural_residual!, x0,
                         method     = :newton,
                         show_trace = true,
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

## Generate DataFrames

# df = DataFrame([ dx θx dy θy dz θz ], :auto)
# rename!(df, [:dx, :θx, :dy, :θy, :dz, :θz])

# data = DataFrame([ Fx dx Mx θx Fy dy My θy Fz dz Mz θz ], :auto)
# rename!(data, [:Fx, :dx, :Mx, :θx, :Fy, :dy, :My, :θy, :Fz, :dz, :Mz, :θz])


## Plotting
#==========================================================================================#

using Plots
gr(dpi = 300)

plot(aspect_ratio = 1)
plot!.(plot_panels(panels[:]), color = :black, label = :none)
plot!(wing_plan, color = :blue, label = :none)
plot!()

hs_pts  = bound_leg_center.(horses)

quiver!(getindex.(hs_pts, 1)[:], getindex.(hs_pts, 2)[:], getindex.(hs_pts, 3)[:], quiver=(getindex.(normies, 1)[:], getindex.(normies, 2)[:], getindex.(normies, 3)[:]), color = :orange)

# new_hs_pts = bound_leg_center.(new_horses)

# plot!.(plot_panels(panels[:]), color = :grey, label = :none)
# quiver!(getindex.(new_hs_pts, 1)[:], getindex.(new_hs_pts, 2)[:], getindex.(new_hs_pts, 3)[:], quiver=(getindex.(new_norms, 1)[:], getindex.(new_norms, 2)[:], getindex.(new_norms, 3)[:]))
# plot!()