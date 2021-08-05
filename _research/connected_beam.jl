using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
# using Einsum

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

# Aerodynamic forces
vlm_forces = surface_forces(aero_surfs[1]) 
sec_forces = sum(vlm_forces, dims = 1)[:] / 2

# Aerodynamic center locations
horsies    = horseshoes(aero_surfs[1])
vlm_acs    = bound_leg_center.(horsies)

# Moments
section_moments(vlm_acs, fem_pts, half_vlm_forces) = sum(permutedims(reduce(hcat, eachrow(vlm_acs) .- Ref(fem_pts))) .× half_vlm_forces, dims = 1)[:]

M_ins  = section_moments(vlm_acs, fem_pts[1:end-1], vlm_forces / 2)
M_outs = section_moments(vlm_acs, fem_pts[2:end],   vlm_forces / 2)

## Load averaging
pt_forces  = adjacent_joiner(sec_forces / 2, sec_forces / 2)
pt_moments = adjacent_joiner(M_ins, M_outs)

## Axis system transformation matrices
#==========================================================================================#

# Compute local beam axes
ss = @. normalize(fem_pts[2:end,:] - fem_pts[1:end-1,:])
ns = @. normalize((plan_coords[end,2:end] - plan_coords[1,1:end-1]) × (plan_coords[1,2:end] - plan_coords[end,1:end-1]))
cs = @. ss × ns

## Active transformations of forces and moments (vs. passive transformation of stiffness matrix)
global_axis     = I(3)            # Global axis system for panels in VLM
beam_local_axis = [ 0. -1.  0. ; 
                    1.  0.  0. ; 
                    0.  0. -1. ]  # Orthogonal local axis system for beam in FEM

local_axes = repeat(beam_local_axis, 1, 1, size(ss)...)

# WTF array of local coordinate systems
wtf = zeros(3, 3, size(ss)...) # ((x,y,z), (c,s,n), chordwise, spanwise)
for inds in CartesianIndices(wtf)
    l,k,i,j  = inds.I
    wtf[l,1,i,j] = cs[i,j][l]
    wtf[l,2,i,j] = ss[i,j][l]
    wtf[l,3,i,j] = ns[i,j][l]
end

## Array comprehension
mats = @views [ local_axes[:,:,inds]' * wtf[:,:,inds]' * global_axis for inds in CartesianIndices(ss) ]

# Einstein summation
# @einsum mats_eins[l,k,i,j] := local_axes[m,l,i,j] * global_axis[m,k]

## Transform forces
fem_forces  = [ mats .* pt_forces[1:end-1] ; [mats[end] * pt_forces[end] ] ]
fem_moments = [ mats .* pt_moments[1:end-1]; [mats[end] * pt_moments[end]] ]
fem_loads   = [ reduce(hcat, fem_forces); reduce(hcat, fem_moments) ]

## Structural setup
#==========================================================================================#

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3 # Density, kg/m³
ν     = 0.3   # Poisson ratio (UNUSED FOR NOW)
R     = 1e-2  # Outer radius, m
t     = 8e-3  # Thickness, m
Ls    = norm.(diff(fem_pts))

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

# Assemble RHS
function assemble_fem_dynamics(Fx, Fy, Fz, Mx, My, Mz)
    px = [ Fx; Mx ]                                 # Fx Mx concatenated
    py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My) # Fy My interspersed
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz) # Fz Mz interspersed

    F  = [ py; pz; px ]
end

function assemble_fem_dynamics(loads)
    px = [ loads[1,:]; loads[4,:] ] # Fx Mx concatenated
    py = loads[[2,5],:][:]          # Fy My interspersed
    pz = loads[[3,6],:][:]          # Fz Mz interspersed
    F  = [ py; pz; px ]
end

Fx = fem_loads[1,:]
Fy = fem_loads[2,:] 
Fz = fem_loads[3,:]
Mx = fem_loads[4,:]
My = fem_loads[5,:]
Mz = fem_loads[6,:]

df_Fs = DataFrame([ Fx Fy Fz Mx My Mz ], :auto)
rename!(df_Fs, [:Fx, :Fy, :Fz, :Mx, :My, :Mz])

## "FEM" setup
K  = tube_stiffness_matrix(aluminum, tubes)

## Loads
Ls = assemble_fem_dynamics(Fx, Fy, Fz, Mx, My, Mz)
# Ls = assemble_fem_dynamics(fem_loads)

## Solve system(s)
xs = K \ Ls

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
df_xs = DataFrame([ dx θx dy θy dz θz ], :auto)
rename!(df_xs, [:dx, :θx, :dy, :θy, :dz, :θz])

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

# Beam circles
n_pts          = 20
circle3D(r, n) = [ (r*cos(θ), 0, r*sin(θ)) for θ in 0:2π/n:2π ];
circ           = circle3D(R, n_pts) 

draw_tube(p1, p2, circ) = [ [ circ_pt .+ p1, circ_pt .+ p2 ] for circ_pt in circ ]

beam_pts   = zip(tupvector(fem_pts[1:end-1]), tupvector(fem_pts[2:end]))
left_pts   = [ draw_tube(pt[1], pt[2], circ) for pt in beam_pts ]

# Beam loads
fem_plot   = reduce(hcat, fem_pts)
loads_plot = [ reduce(hcat, pt_forces); reduce(hcat, pt_moments) ]

# Aerodynamic centers and forces
panel_plot = plot_panels(panels[:])
ac_plot    = reduce(hcat, vlm_acs)
force_plot = reduce(hcat, vlm_forces)

# Beam axes
axes_plot  = reduce(hcat, (fem_pts[1:end-1] + fem_pts[2:end]) / 2)
cs_plot    = reduce(hcat, cs)
ss_plot    = reduce(hcat, ss)
ns_plot    = reduce(hcat, ns)

# Plot
b = aero_state.span_ref
plot(camera = (45, 45), 
     xlim = (-b/2, b/2),
    #  ylim = (-b/2, b/2), 
     zlim = (-b/2, b/2)
    )

# Panels
[ plot!(pans, color = :black, label = ifelse(i == 1, "Panels", :none)) for (i, pans) in enumerate(panel_plot) ]

# Planform
plot!(wing_plan, color = :blue, label = "Planform")

# Beams
[ plot!(reduce(vcat, pt), color = :green, label = ifelse(i == 1, "Beams", :none)) for (i, pt) in enumerate(left_pts) ]

# Forces
quiver!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:], quiver=(loads_plot[1,:], loads_plot[2,:], loads_plot[3,:] ) .* 0.1, label = "Beam Forces")
quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1, label = "Panel Forces")
scatter!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], label = "Aerodynamic Centers")

# Axis systems
quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:], quiver=(cs_plot[1,:], cs_plot[2,:], cs_plot[3,:]) .* 1e-1, color = :orange, label = :none)
quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:], quiver=(ns_plot[1,:], ns_plot[2,:], ns_plot[3,:]) .* 1e-1, color = :red,    label = :none)
quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:], quiver=(ss_plot[1,:], ss_plot[2,:], ss_plot[3,:]) .* 1e-1, color = :brown,  label = :none)

plot!()