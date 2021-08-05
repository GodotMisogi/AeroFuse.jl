using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames

## Aerodynamic setup
#==========================================================================================#

# Define wing
wing = WingSection(root_foil  = naca4(0,0,1,5),
                   span       = 3.11,
                   dihedral   = 0.0, # 5.0
                   sweep_LE   = 0.0, # 20.0
                   taper      = 1.0, # 0.5
                   root_chord = 0.3, 
                   root_twist = 0.0, # 2.0
                   tip_twist  = 0.0)
wing_mac    = mean_aerodynamic_center(wing)
wing_plan   = plot_wing(wing)
wing_name   = "Wing"
print_info(wing, wing_name)

# Mesh
span_num        = 5
chord_num       = 2
panels, normies = panel_wing(wing, span_num, chord_num, spacing = "sine");
aircraft        = Dict(wing_name => (panels, normies));

# Set up aerodynamic state
aero_state = VLMState(1., 0., 0., [0.0, 0.0, 0.0], 
                      rho_ref   = 1.225,
                      r_ref     = [ wing_mac[1], 0., 0. ],
                      area_ref  = projected_area(wing), 
                      chord_ref = mean_aerodynamic_chord(wing), 
                      span_ref  = span(wing));

# Test case - Fixed speed
aero_state.speed   = 22.876
aero_state.alpha   = deg2rad(5.)
aero_state.beta    = deg2rad(0.)
aero_state.rho_ref = 0.770816

# Build system with initial guess from aerodynamic-only analysis
aero_system = solve_case(aircraft, aero_state)
aero_surfs  = surfaces(aero_system)
print_coefficients(aero_surfs[1], aero_state);

## Load transfer scheme
#==========================================================================================#

# Functions on adjacencies
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

## Mesh setup
fem_w       = 0.35
plan_coords = chord_coordinates(wing, [span_num], chord_num, spacings = ["sine"])
mesh        = reshape(permutedims(reduce(hcat, plan_coords)), (size(plan_coords)..., 3))

# FEM beam node locations as Matrix
fem_pts  = (1 - fem_w) * plan_coords[1,:] + fem_w * plan_coords[end,:]

# Aerodynamic center locations and forces
forces      = surface_forces(aero_surfs[1]) 
half_forces = sum(forces, dims = 1)[:] / 2
horsies     = horseshoes(aero_surfs[1])
ac_mine     = bound_leg_center.(horsies)

# Moments
M_ins  = sum(permutedims(reduce(hcat, eachrow(ac_mine) .- Ref(fem_pts[1:end-1]))) .× forces / 2, dims = 1)[:]
M_outs = sum(permutedims(reduce(hcat, eachrow(ac_mine) .- Ref(fem_pts[2:end])))   .× forces / 2, dims = 1)[:]

## Load averaging
pt_forces   = adjacent_joiner(half_forces, half_forces)
pt_moments  = adjacent_joiner(M_ins, M_outs)

## Axis system transformation matrices
#==========================================================================================#

# Compute local beam axes
ss = @. normalize(fem_pts[2:end,:] - fem_pts[1:end-1,:])
ns = @. normalize((plan_coords[end,2:end] - plan_coords[1,1:end-1]) × (plan_coords[1,2:end] - plan_coords[end,1:end-1]))
cs = @. ss × ns

# WTF array version
wtf = zeros(3, 3, size(ss)...) # ((x,y,z), (c,s,n), chordwise, spanwise)
for inds in CartesianIndices(wtf)
    l,k,i,j  = inds.I
    wtf[l,1,i,j] = cs[i,j][l]
    wtf[l,2,i,j] = ss[i,j][l]
    wtf[l,3,i,j] = ns[i,j][l]
end

## Active transformations of forces and moments (vs. passive transformation of stiffness matrix)
global_axis     = I(3)        # Global axis system for panels in VLM
beam_local_axis = [0. 1. 0.; 
                   0. 0. 1.; 
                   1. 0. 0.]  # Orthogonal local axis system for beam in FEM

local_ref = repeat(beam_local_axis, 1, 1, size(ss)...)

## Array comprehension
mats = [ local_ref[:,:,inds]' * wtf[:,:,inds] * global_axis for inds in CartesianIndices(ss) ]

# Einstein summation
# @einsum mats_eins[l,k,i,j] := local_ref[m,l,i,j] * wtf[m,n,i,j] * global_axis[n,k]

## Transform forces
force_trans = mats .* forces

## Matrix version
# ac_fem   = reshape(permutedims(reduce(hcat, ac_mine)), (size(ac_mine)..., 3))
# fem_arr  = (1 - fem_w) * mesh[1,:,:] + fem_w * mesh[end,:,:]

# Moments
# M_in  = [ let xs = (ac_fem[i,:,:] - fem_arr[1:end-1,:]); 
#               @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
#               for i in eachindex(ac_fem[:,1,1]) ]
# M_out = [ let xs = (ac_fem[i,:,:] - fem_arr[2:end,:]); 
#               @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
#               for i in eachindex(ac_fem[:,1,1]) ]

# M_ins  = sum(M_in)
# M_outs = sum(M_out) 


fem_loads  = zeros(6, length(pt_forces))
fem_forces = reduce(hcat, half_forces)
fem_loads[1:3,1:end-1] = fem_forces
fem_loads[1:3,2:end]  += fem_forces

fem_M_ins  = reduce(hcat, M_ins)
fem_M_outs = reduce(hcat, M_outs)
fem_loads[4:end,1:end-1] = fem_M_ins
fem_loads[4:end,2:end]  += fem_M_outs

# Constrain the plane of symmetry and print
# fem_loads[:,span_num+1] .= 0.
fem_loads
 
## Transform loads to principal axes
fem_force_trans  = [ dircos[1] * fem_loads[1:3,1:span_num]   dircos[2] * fem_loads[1:3,span_num+1:end]   ]
fem_moment_trans = [ dircos[1] * fem_loads[4:end,1:span_num] dircos[2] * fem_loads[4:end,span_num+1:end] ]
fem_loads_trans  = [ fem_force_trans  ; 
                     fem_moment_trans ]
fem_loads_trans

## Splitting for Wing case
#==========================================================================================#

function middle_index(x :: AbstractArray)
    n = length(x)
    if n % 2 == 0
        Int(n / 2)
    else
        ceil(Int, n / 2)
    end
end

middle(x :: AbstractVector) = @view x[middle_index(x)]

zero_vec = [SVector(0,0,0.)]

# Ls    = (norm ∘ bound_leg_vector).(horsies)
# Fs  = [ [ zero_vec; pt_forces[2:end]  ] ]
# Ms  = [ [ zero_vec; pt_moments[2:end] ] ]
# F_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Fs))
# M_S = map(x -> [ x[1] * vec for vec in x[2] ], zip(dircos, Ms))

Ls    = norm.(diff(fem_pts))

# NEED TO CHECK THE DIRECTIONS
left_forces   = @views pt_forces[1:middle_index(pt_forces)-1]
left_moments  = @views pt_moments[1:middle_index(pt_moments)-1]
right_forces  = @views [ zero_vec; pt_forces[middle_index(pt_forces)+1:end] ]
right_moments = @views [ zero_vec; pt_moments[middle_index(pt_moments)+1:end] ]


## Axis transformation of loads
F_S = map(x -> map(y -> x[1] * y, x[2]), zip(dircos, [left_forces,  right_forces ]))
M_S = map(x -> map(y -> x[1] * y, x[2]), zip(dircos, [left_moments, right_moments]))

## Structural setup
#==========================================================================================#

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3   # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-2  # Outer radius, m
t     = 8e-3  # Thickness, m

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

# Assemble RHS
function assemble_fem_dynamics(pt_forces, pt_moments)
    Fx = getindex.(pt_forces,  1)
    Fy = getindex.(pt_forces,  2) 
    Fz = getindex.(pt_forces,  3)
    Mx = getindex.(pt_moments, 1)
    My = getindex.(pt_moments, 2)
    Mz = getindex.(pt_moments, 3)

    px = (collect ∘ Iterators.flatten ∘ zip)(Fx, Mx)
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
    py = [ Fy[:]; My[:] ]

    F  = [ pz; px; py ]
end

function assemble_fem_dynamics(loads)
    px = loads[[2,5],:][:]          # Fy My interspersed
    pz = loads[[3,6],:][:]          # Fz Mz interspersed
    py = [ loads[1,:]; loads[4,:] ] # Fx Mx concatenated
    F  = [ pz; px; py ]
end

Fx = getindex.(pt_forces,  1)
Fy = getindex.(pt_forces,  2) 
Fz = getindex.(pt_forces,  3)
Mx = getindex.(pt_moments, 1)
My = getindex.(pt_moments, 2)
Mz = getindex.(pt_moments, 3)

df = DataFrame([ Fx Fy Fz Mx My Mz ], :auto)
rename!(df, [:Fx, :Fy, :Fz, :Mx, :My, :Mz])

# py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My)
# pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz)
# px = [ Fx[:]; Mx[:] ]

## "FEM" setup
K  = tube_stiffness_matrix(aluminum, tubes)

## Loads
Fs = assemble_fem_dynamics.(F_S, M_S)
# F  = reduce(vcat, Fs)
F  = assemble_fem_dynamics(fem_loads_trans) 

## Solve system(s)
xs = K \ F

## Displacement transfer scheme
#==========================================================================================#

function compute_displacements(Δs, n)
    # Assembly
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

    dx, θx, dy, θy, dz, θz
end

dx, θx, dy, θy, dz, θz = compute_displacements(xs, length(fem_pts))
ds, θs = SVector.(dx, dy, dz), SVector.(θx, θy, θz)

## Generate DataFrames
df = DataFrame([ dx θx dy θy dz θz ], :auto)
rename!(df, [:dx, :θx, :dy, :θy, :dz, :θz])

## Testing permutations and transformations of stiffness matrices
#==========================================================================================#

# Testing blocks
Ks = [ tube_stiffness_matrix(aluminum, [tube]) for tube in tubes ]
Ks[1]

## Testing simultaneous permutations of rows and columns:
# 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
# 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]
inds = [9,1,5,11,6,2,10,3,7,12,8,4]
perm_Ks = [ K[inds,:][:,inds] for K in Ks ] 
perm_Ks[1]

## Axis transformation of stiffness matrices
axis_trans = [0 1 0; 0 0 -1; -1 0 0]
K_trans = kron(I(4), axis_trans)

K_tran = repeat(K_trans, 1, 1, length(Ls))

perm_K = zeros(12, 12, length(Ls))
for i in 1:length(Ls)
    perm_K[:,:,i] .= perm_Ks[i]
end

## Using Einstein summation convention for multiplication: D = TᵗKT
using Einsum
@einsum D[j,k,i] := K_tran[l,j,i] * perm_K[l,m,i] * K_tran[m,k,i]

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

n_pts = 20
circle3D(r, n) = [ (r*cos(θ), 0, r*sin(θ)) for θ in 0:2π/n:2π ];
circ     = circle3D(R, n_pts) 

beam_pts = collect(zip(tupvector(fem_pts[1:end-1]), tupvector(fem_pts[2:end])))
left_pts = [ [ [ circ_pt .+ pt[1], circ_pt .+ pt[2] ] for circ_pt in circ ] for pt in beam_pts ]

# Aerodynamic forces
mid_pans = ac_mine # midpoint.(panels)[:]
mid_xs   = getindex.(mid_pans, 1)
mid_ys   = getindex.(mid_pans, 2)
mid_zs   = getindex.(mid_pans, 3)

# Beam midpoints
fem_mids = (fem_pts[1:end-1] + fem_pts[2:end]) / 2
fem_xs   = getindex.(fem_mids, 1)
fem_ys   = getindex.(fem_mids, 2)
fem_zs   = getindex.(fem_mids, 3)

# Forces between panels
cs_plot    = reduce(hcat, cs)
ss_plot    = reduce(hcat, ss)
ns_plot    = reduce(hcat, ns)
force_plot = reduce(hcat, forces)

#
hs_pts = reduce(hcat, points.(eachrow(bound_leg.(horsies))))
hs_xs  = getindex.(hs_pts, 1)
hs_ys  = getindex.(hs_pts, 2)
hs_zs  = getindex.(hs_pts, 3)

# Plot
b = aero_state.span_ref
plot(camera = (45, 45), 
     xlim = (-b/2, b/2),
    #  ylim = (-b/2, b/2), 
     zlim = (-b/2, b/2)
    )

# Panels
[ plot!(pans, color = :black, label = ifelse(i == 1, "Panels", :none)) for (i, pans) in enumerate(plot_panels(panels[:])) ]

# Planform
plot!(wing_plan, color = :blue, label = "Planform")

# Beams
[ plot!(reduce(vcat, pt), color = :green, label = ifelse(i == 1, "Beams", :none)) for (i, pt) in enumerate(left_pts) ]

# Forces
quiver!(fem_arr[:,1], fem_arr[:,2], fem_arr[:,3], quiver=(fem_loads[1,:],  fem_loads[2,:],  fem_loads[3,:] ) .* 0.1, label = "Beam Forces")
quiver!(mid_xs[:], mid_ys[:], mid_zs[:], quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1, label = "Panel Forces")

# Axis systems
quiver!(fem_xs[:], fem_ys[:], fem_zs[:], quiver=(cs_plot[1,:], cs_plot[2,:], cs_plot[3,:]) .* 1e-1, color = :orange, label = :none)
quiver!(fem_xs[:], fem_ys[:], fem_zs[:], quiver=(ns_plot[1,:], ns_plot[2,:], ns_plot[3,:]) .* 1e-1, color = :red,    label = :none)
quiver!(fem_xs[:], fem_ys[:], fem_zs[:], quiver=(ss_plot[1,:], ss_plot[2,:], ss_plot[3,:]) .* 1e-1, color = :brown,  label = :none)

# Nodes
scatter!(ac_fem[:,:,1][:], ac_fem[:,:,2][:], ac_fem[:,:,3][:], label = "Aerodynamic Centers")
# scatter!(getindex.(ac_mine, 1)[:], getindex.(ac_mine, 2)[:], getindex.(ac_mine, 3)[:], label = "ACs (Mine)")
plot!()