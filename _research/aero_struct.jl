##
using Revise
using AeroMDAO
using Base.Iterators
using LinearAlgebra
using StaticArrays
# using ForwardDiff
using DataFrames
using NLsolve

## Helper functions
#==========================================================================================#

# Sum adjacent values
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

## Aerodynamic analysis
#==========================================================================================#

aerodynamic_forces(surfs) = reduce(vcat, (vec ∘ surface_forces).(surfs))

points(horses) = @. r1(bound_leg_center(horses), horses)[:], r2(bound_leg_center(horses), horses)[:]

lift_force(surfs) = sum(sum ∘ surface_forces, surfs) 

function solve_aerodynamic_residual!(R, Γ, system :: VLMSystem, state :: VLMState)
    V = freestream_to_cartesian(-state.U, state.alpha, state.beta)

    # Assemble matrix system
    compute_influence_matrix!(system, V)
    compute_boundary_condition!(system, V, state.omega)

    # Evaluate residual
    evaluate_residual!(R, Γ, system)

    R
end

function update_circulations!(Γ, surfs :: Vector{<: VLMSurface}) 
    # Get sizes and indices for reshaping (UGLY AF)
    sizes = @. (size ∘ horseshoes)(surfs)
    inds  = [ 0; cumsum(prod.(sizes)) ]
    Γs    = reshape_array(Γ, inds, sizes)

    # Allocate surface circulations
    for (surf, Γ_vec) in zip(surfs, Γs)
        surf.circulations = Γ_vec
    end

    surfs
end

function compute_loads(forces, r1s, r2s)
    # Aerodynamic forces
    half_forces = forces ./ 2
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

solve_beam_residual!(R, K, δ, F) = 
    R .= K * δ - F

function compute_displacements(Δs, r1s, r2s)
    # Assembly
    n  = length(horses) + 1 # Number of finite element points
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

function transfer_displacements!(system, xyzs, normals, ds, θs)
    # Displace and rotate about quarter-chord point, where the beam is located. Needs many corrections...
    Δxyzs      = permutedims([ ds ds ])
    new_xyzs   = xyzs + Δxyzs
    new_pans   = make_panels(new_xyzs)
    new_horses = horseshoe_line.(new_pans) 
    new_norms  = @. normals[:] + (θs[1:end-1] + θs[2:end]) / 2 # WRONG

    # Set new system variables
    system.horseshoes = new_horses[:]
    system.normals    = new_norms[:]

    new_xyzs
end

## Weights and fuel loads
#==========================================================================================#

evaluate_weight_residual(L, W, n) = L - n * W

## Coupled residual system
#==========================================================================================#

function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_state :: VLMState, surfs, xyzs, normals, K, W, load_factor)
    n = (prod ∘ size)(horseshoes(system))   # Get panel size

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]
    α = @view x[end]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute and transfer displacements
    r1s, r2s = (points ∘ horseshoes)(system)
    ds, θs = compute_displacements(δ, r1s, r2s)
    transfer_displacements!(system, xyzs, normals, ds, θs)

    # Aerodynamic residuals
    aero_state.alpha = α[1]
    solve_aerodynamic_residual!(R_A, Γ, aero_system, aero_state)
    update_circulations!(Γ, surfs)
    
    # Compute forces
    forces = aerodynamic_forces(surfs)
    F      = compute_loads(forces, r1s, r2s)
    L      = lift_force(surfs)[3]

    # Structural residuals
    solve_beam_residual!(R_S, K, δ, F)

    # Weight residual
    R_W = evaluate_weight_residual(L, W, load_factor)

    R
end

## Case
#==========================================================================================#

## Aerodynamic variables

# Define wing
wing = WingSection(span       = 4.0,
                   dihedral   = 0.0,
                   sweep_LE   = 15.0,
                   taper      = 0.4,
                   root_chord = 2.0,
                   root_twist = 0.0,
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
name        = "Wing"
print_info(wing)

# Mesh
span_num        = 5
chord_num       = 1
xyzs            = chord_coordinates(wing, [span_num], chord_num; spacings = ["cosine"])
panels, normals = panel_wing(wing, span_num, chord_num);
aircraft        = Dict(name => (panels, normals));

## Define VLMSystem

# Set up state
state = VLMState(Freestream(1.0, 1.0, 0.0, [0.0, 0.0, 0.0]), 
                 rho_ref   = 1.225,
                 r_ref     = SVector(wing_mac[1], 0., 0.),
                 area_ref  = projected_area(wing), 
                 chord_ref = mean_aerodynamic_chord(wing), 
                 span_ref  = span(wing));

# Build system
system, surfs = build_system(aircraft);
horses  = horseshoes(system)
normies = system.normals 

## Weight variables

# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]
weight      = 10.
load_factor = 2

## Structural variables

# Tube properties
E  = 50e9
G  = 30e9
J  = 2.
A  = 0.1
Iy = 5
Iz = 5
Ls = @. (norm ∘ bound_leg_vector)(horses)[:]
n  = length(Ls)
x  = [ fill(E, n) fill(G, n) fill(A, n) fill(Iy, n) fill(Iz, n) fill(J, n) Ls ]

# Stiffness matrix setup
K = tube_stiffness_matrix(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7])

## Coupled system evaluation
vec_size = (prod ∘ size ∘ horseshoes)(system) + (n + 1) * 6 + 1

R   = zeros(vec_size)
x   = rand(vec_size)
eva = solve_coupled_residual!(R, x, system, state, surfs, xyzs, normies, K, weight, load_factor)

## Coupled system solution
solve_aerostructural_residual!(R, x) = solve_coupled_residual!(R, x, system, state, surfs, xyzs, normies, K, weight, load_factor)

res = nlsolve(solve_aerostructural_residual!, x,
            #   method   = :newton,
            #   autodiff = :forward,
              show_trace = true,
              )

##
# df = DataFrame([ dx θx dy θy dz θz ], :auto)
# rename!(df, [:dx, :θx, :dy, :θy, :dz, :θz])

data = DataFrame([ Fx dx Mx θx Fy dy My θy Fz dz Mz θz ], :auto)
rename!(data, [:Fx, :dx, :Mx, :θx, :Fy, :dy, :My, :θy, :Fz, :dz, :Mz, :θz])


## Plotting
#==========================================================================================#

using Plots
# unicodeplots()
# spy(K)
hs_pts     = bound_leg_center.(horses)
new_hs_pts = bound_leg_center.(new_horses)


gr(dpi = 300)
plot()
plot!.(plot_panels(new_pans[:]), color = :black, label = :none)

quiver!(getindex.(hs_pts, 1)[:], getindex.(hs_pts, 2)[:], getindex.(hs_pts, 3)[:], quiver=(getindex.(normies, 1)[:], getindex.(normies, 2)[:], getindex.(normies, 3)[:]), color = :orange)

plot!.(plot_panels(panels[:]), color = :grey, label = :none)
quiver!(getindex.(new_hs_pts, 1)[:], getindex.(new_hs_pts, 2)[:], getindex.(new_hs_pts, 3)[:], quiver=(getindex.(new_norms, 1)[:], getindex.(new_norms, 2)[:], getindex.(new_norms, 3)[:]))
plot!()