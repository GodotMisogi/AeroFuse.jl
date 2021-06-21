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
begin
    state = VLMState(fs, 
                     rho_ref   = ρ,
                     r_ref     = SVector(ref...),
                     area_ref  = S, 
                     chord_ref = c, 
                     span_ref  = b);

    system, surfs, nf_t, ff_t = solve_case!(aircraft, state);
    print_coefficients(nf_t, ff_t, state.name)
end

## Evaluate residual
solve_aerodynamics!(R, x) = solve_residual!(system, R, x)

Γ_0 = rand(span_num * chord_num)
res = nlsolve(solve_aerodynamics!, Γ_0, 
              show_trace = true,
              autodiff   = :forward) # Works correctly, unsurprisingly

## Processing for structures
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

q           = dynamic_pressure(ρ, V)
areas       = @. panel_area(horseshoe_panels)[:]
half_forces = @. force(CFs[:], q, S) / 2
r1s         = @. r1(bound_leg_center(horseshoes), horseshoes)[:]
r2s         = @. r2(bound_leg_center(horseshoes), horseshoes)[:]
M1s         = @. r1s × half_forces
M2s         = @. r2s × half_forces
zero_vec    = SVector(0., 0., 0.)
weight_vec  = SVector(0., 0., 1.)

# Interspersing
forces            = adjacent_joiner(half_forces, half_forces)

n                 = ceil(Int, length(horseshoe_panels) / 2) + 1
forces[n-1:n+1] .-= forces[n-1:n+1]   # Boundary condition, setting F = 0 at center of wing

moments           = adjacent_joiner(M1s, M2s)

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

## Tube properties
E  = 50e9
G  = 30e9
J  = 2.
A  = 0.1
Iy = 5.
Iz = 5.
Ls = @. (norm ∘ bound_leg_vector)(horseshoes)[:]
n  = length(Ls)
x  = [ fill(E, n) fill(G, n) fill(A, n) fill(Iy, n) fill(Iz, n) fill(J, n) Ls ]

## Stiffness matrix setup
K = tube_stiffness_matrix(x[:,1], x[:,2], x[:,3], x[:,4], x[:,5], x[:,6], x[:,7])

## Solve system
f = ps
δ = K \ f

## Get displacements
n1 = length(Fy) + length(My)
n2 = n1 + length(Fz) + length(Mz)
n3 = n2 + length(Fx)
n4 = n3 + length(Mx)

dy = δ[1:2:n1]
θy = δ[2:2:n1]
dz = δ[n1+1:2:n2]
θz = δ[n1+2:2:n2]
dx = δ[n2+1:n3]
θx = δ[n3+1:n4]

state = DataFrame([ dx θx dy θy dz θz ], :auto)
rename!(state, [:dx, :θx, :dy, :θy, :dz, :θz])

## Transfer displacements to aerodynamic mesh

# Displace and rotate about quarter-chord point, where the beam is located. Needs many corrections...
rs        = @. SVector(dx, dy, dz)
θs        = @. SVector(θx, θy, θz)
ds        = rs .+ θs .× adjacent_joiner(r1s, r2s)
xyzs      = chord_coordinates(wing, [Int(span_num / 2)], chord_num; spacings = ["cosine"])
new_xyzs  = xyzs + permutedims([ ds ds ])
new_pans  = make_panels(new_xyzs)
new_norms = @. normals[:] + (θs[1:end-1] + θs[2:end]) / 2 # WRONG

# A bijective mapping between the wing's geometric and coordinate representations needs to be defined, if it exists.
# Suspecting it doesn't due to different perturbations in the section spans at the leading and trailing edges.

##
data = DataFrame([ Fx dx Mx θx Fy dy My θy Fz dz Mz θz ], :auto)
rename!(data, [:Fx, :dx, :Mx, :θx, :Fy, :dy, :My, :θy, :Fz, :dz, :Mz, :θz])

## Plotting
using Plots
# unicodeplots()
# spy(K)

plotly(dpi = 300)
plot(aspect_ratio = 1)
plot!.(plot_panels(new_pans[:]), color = :black, label = "new")
plot!.(plot_panels(horseshoe_panels[:]), color = :grey, label = "old")
plot!()  