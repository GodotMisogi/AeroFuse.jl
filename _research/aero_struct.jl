##
using Revise
using AeroMDAO
using Base.Iterators
using LinearAlgebra
using StaticArrays
using ForwardDiff
using DataFrames

## Define wing
wing = WingSection(span = 4.0, 
                   dihedral = 1.0, 
                   sweep_LE = 15.0, 
                   taper = 0.4, 
                   root_chord = 2.0, 
                   root_twist = 0.0, 
                   tip_twist = 0.0)
wing_mac  = mean_aerodynamic_center(wing)
wing_plan = plot_wing(wing)
print_info(wing)

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

@time nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = 
    solve_case(wing, fs; 
               rho_ref   = ρ, 
               r_ref     = ref,
               span_num  = span_num, 
               chord_num = chord_num,
               viscous   = true,
               x_tr      = 0.3);
               
print_coefficients(nf_coeffs, ff_coeffs, "Wing")

## Needs to be moved to NonDimensional
force(CF, q, S) = CF * q * S
moment(CM, q, S, c) = CM * q * S * c
moment(CM, q, S, b, c) = moment.(CM, q, S, [b, c, b])

## Processing for structures
q           = dynamic_pressure(ρ, V)
areas       = @. panel_area(horseshoe_panels)[:]
half_forces = @. force(CFs[:], q, S) / 2
M1s         = @. r1(bound_leg_center(horseshoes), horseshoes)[:] × half_forces
M2s         = @. r2(bound_leg_center(horseshoes), horseshoes)[:] × half_forces 

weight_vector = SVector(0., 0., 1.)

# Interspersing

forces            = [ half_forces; [SVector(0.,0.,0.)] ] .+ [ [SVector(0.,0.,0.)]; half_forces ]

n                 = ceil(Int, length(horseshoe_panels) / 2) + 1
forces[n-1:n+1] .-= forces[n-1:n+1]   # Boundary condition, setting F = 0 at center of wing

moments           = [         M1s; [SVector(0.,0.,0.)] ] .+ [ [SVector(0.,0.,0.)]; M2s ]

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

##
data = DataFrame([ Fx dx Mx θx Fy dy My θy Fz dz Mz θz ], :auto)
rename!(data, [:Fx, :dx, :Mx, :θx, :Fy, :dy, :My, :θy, :Fz, :dz, :Mz, :θz])

## Plotting
using Plots
unicodeplots()
spy(K)

# plotly(dpi = 300)
# plot(wing_plan, aspect_ratio = 1, zlims = (-0.1, 1.0))