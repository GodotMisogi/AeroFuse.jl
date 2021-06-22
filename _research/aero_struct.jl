##
using Revise
using AeroMDAO
using Base.Iterators
using LinearAlgebra
using StaticArrays
# using ForwardDiff
using DataFrames
using NLsolve

## Define and mesh wing
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
panels, normals = panel_wing(wing, 5, 1);
aircraft = Dict(name => (panels, normals));

## Define freestream variables and reference values
ρ 		= 1.225
ref 	= [wing_mac[1], 0., 0.]
V, α, β = 1.0, 1.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 		= Freestream(V, α, β, Ω)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

## Evaluate case
span_num  = 10
chord_num = 1

# Set up state
state = VLMState(fs, 
                 rho_ref   = ρ,
                 r_ref     = SVector(ref...),
                 area_ref  = S, 
                 chord_ref = c, 
                 span_ref  = b);

# Solve system
system, surfs, coeffs = solve_case!(aircraft, state);
nf_t, ff_t = coeffs[state.name]
print_coefficients(nf_t, ff_t, state.name)

# Get relevant variables
CFs    = surface_force_coefficients(surfs[name], state)
horses = horseshoes(surfs[name])

normies = system.normals # Need to check why the getter function isn't working

## Alternative: Evaluate residual
solve_aerodynamics!(R, x) = solve_residual!(system, R, x)

Γ_0 = rand(span_num * chord_num)
res = nlsolve(solve_aerodynamics!, Γ_0, 
              show_trace = true,
              autodiff   = :forward) # Works correctly, unsurprisingly

## Processing for structures
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

q           = dynamic_pressure(ρ, V)
areas       = @. panel_area(panels)[:]
forces      = @. force(CFs[:], q, S)
r1s         = @. r1(bound_leg_center(horses), horses)[:]
r2s         = @. r2(bound_leg_center(horses), horses)[:]

function load_transfer(forces, r1s, r2s)
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

    ## Assembly
    Fx = getindex.(forces, 1)
    Fy = getindex.(forces, 2)
    Fz = getindex.(forces, 3)
    Mx = getindex.(moments, 1)
    My = getindex.(moments, 2)
    Mz = getindex.(moments, 3)

    py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
    px = [ Fx[:]; Mx[:] ]

    ps = [ py; pz; px ]
    F_A = ps
end

F_A = load_transfer(forces, r1s, r2s)

## Weights and fuel loads
# W   = force(ff_t[3], q, S)
# W_y = fill(W / length(CFs) + 1, length(CFs) + 1)
# W   = (collect ∘ Iterators.flatten ∘ zip)(W_y, zeros(length(My)))
# F_W = [ zeros(length(py)); W; zeros(length(px)) ]

solve_weight_residual!(R_W, L, W, n) = R_W .= L - n * W

## Tube properties
E  = 50e9
G  = 30e9
J  = 2.
A  = 0.1
Iy = 5
Iz = 5
Ls = @. (norm ∘ bound_leg_vector)(horses)[:]
n  = length(Ls)
x  = [ fill(E, n) fill(G, n) fill(A, n) fill(Iy, n) fill(Iz, n) fill(J, n) Ls ]

## Stiffness matrix setup
K = tube_stiffness_matrix(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7])

## Solve system
F = F_A #- F_W
δ = K \ F

## Get displacements
n  = length(horses) + 1 # Number of finite element points
n1 = 2n                 # dim(Fy + My)
n2 = n1 + 2n            # dim(Fz + Mz)
n3 = n2 +  n            # dim(Fx)
n4 = n3 +  n            # dim(Mx)

dy = δ[1:2:n1]
θy = δ[2:2:n1]
dz = δ[n1+1:2:n2]
θz = δ[n1+2:2:n2]
dx = δ[n2+1:n3]
θx = δ[n3+1:n4]

df = DataFrame([ dx θx dy θy dz θz ], :auto)
rename!(df, [:dx, :θx, :dy, :θy, :dz, :θz])

## Transfer displacements to aerodynamic mesh

# Displace and rotate about quarter-chord point, where the beam is located. Needs many corrections...
rs        = @. SVector(dx, dy, dz)
θs        = @. SVector(θx, θy, θz)
ds        = rs .+ θs .× adjacent_joiner(r1s, r2s)
xyzs      = chord_coordinates(wing, [Int(span_num / 2)], chord_num; spacings = ["cosine"])
δs        = permutedims([ ds ds ])
new_xyzs  = xyzs + δs
new_pans  = make_panels(new_xyzs)
new_horses = horseshoe_line.(new_pans) 
new_norms = @. normals[:] + (θs[1:end-1] + θs[2:end]) / 2 # WRONG

# Set new normals
system.normals = new_norms

# Solve new system with displaced coordinates and normals, obviously wrong for now
solve_aerodynamics(Γs) = solve_residual(system, state, xyzs, δs, Γs)
res = nlsolve(solve_aerodynamics, Γ_0)

# Allocate new circulations
system.circulations = res.zero

# A bijective mapping between the wing's geometric and coordinate representations needs to be defined, if it exists.
# Suspecting it doesn't due to different perturbations in the section spans at the leading and trailing edges.

## Coupled residuals
function solve_coupled_residual!(R, x, aero_state, aircraft, K, F, W, load_factor)
    n = (prod ∘ size)(aircraft["Wing"][1])   # Get panel size

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]
    α = x[end]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = [ R[end] ]

    # Solve system and get forces
    aero_state.alpha = α
    solve_residual!(aero_system, R_A, Γ)

    # aero_system, surfs, coeffs = solve_case!(aircraft, aero_state);

    # (Wrong setup because this gets the aero-only force, needs rethinking)
    L = force(coeffs[aero_state.name][2][3], dynamic_pressure(aero_state.rho_ref, aero_state.U), aero_state.area_ref)

    # Aerodynamic residuals

    # Structural residuals
    AeroMDAO.Beams.solve_beam_residual!(R_S, K, δ, F)

    # Weight residual
    solve_weight_residual!(R_W, L, W, load_factor)

    return R
end

vec_size = (span_num * chord_num) + (span_num + 1) * 6 + 1

R   = zeros(vec_size)
x   = rand(vec_size)
res = solve_coupled_residual!(R, x, state, aircraft, K, F_A, 10., 2.)

##
data = DataFrame([ Fx dx Mx θx Fy dy My θy Fz dz Mz θz ], :auto)
rename!(data, [:Fx, :dx, :Mx, :θx, :Fy, :dy, :My, :θy, :Fz, :dz, :Mz, :θz])

## Plotting
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