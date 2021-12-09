using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays
using DataFrames
using Einsum
using UnicodePlots

## Aerodynamic setup
#==========================================================================================#

# Define wing
wing = WingSection(root_foil  = naca4(0,0,1,5),
                   span       = 3.11,
                   dihedral   = 0.0, # 5.0
                   sweep_LE   = 20.0, # 20.0
                   taper      = 0.5, # 0.5
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
panels, normies = panel_wing(wing, span_num, chord_num, spacing = Cosine());
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

## Mesh setup
vlm_mesh = chord_coordinates(wing, [span_num], chord_num, spacings = [Cosine()])

## Aerodynamic forces
vlm_forces = surface_forces(aero_surfs[1])
sec_forces = sum(vlm_forces, dims = 1)[:] / 2

## Aerodynamic center locations
horsies    = horseshoes(aero_surfs[1])
vlm_acs    = bound_leg_center.(horsies)

## This is lazy?
section_moments(vlm_acs, fem_mesh, half_vlm_forces) = sum(x -> (x[1] .- fem_mesh) .× x[2], zip(eachrow(vlm_acs), eachrow(half_vlm_forces)))

## Load averaging
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

## FEM beam node locations as Matrix
make_beam_mesh(vlm_mesh, fem_w = 0.35) =  (1 - fem_w) * vlm_mesh[1,:] + fem_w * vlm_mesh[end,:]

function compute_loads(vlm_acs, vlm_forces, fem_mesh)
    # Forces
    sec_forces   = sum(vlm_forces, dims = 1)[:] / 2
    beam_forces  = adjacent_joiner(sec_forces / 2, sec_forces / 2)

    # Moments
    M_ins        = @views section_moments(vlm_acs, fem_mesh[1:end-1], vlm_forces / 2)
    M_outs       = @views section_moments(vlm_acs, fem_mesh[2:end],   vlm_forces / 2)
    beam_moments = adjacent_joiner(M_ins, M_outs)

    # Concatenate forces and moments into loads array
    fem_loads    = [ reduce(hcat, beam_forces); reduce(hcat, beam_moments) ]
end

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

# Make FEM mesh
fem_w     = 0.35
fem_mesh  = make_beam_mesh(vlm_mesh, fem_w)
Ls        = norm.(diff(fem_mesh))

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

## Axis transformation of stiffness matrices
function axis_transformation(fem_mesh)
    # Compute local beam axes
    ss = @. normalize(fem_mesh[2:end] - fem_mesh[1:end-1])
    ns = ss .× fill(SVector(1,0,0), length(fem_mesh) - 1)
    cs = @. ss × ns

    # WTF array of local coordinate systems - I can do better than this -_-
    wtf = zeros(3, 3, size(ss)...) # ((x,y,z), (c,s,n), chordwise, spanwise)
    for inds in CartesianIndices(wtf)
        l,k,i  = inds.I
        wtf[l,1,i] =  ns[i][l]
        wtf[l,2,i] = -cs[i][l]
        wtf[l,3,i] = -ss[i][l]
    end
    wtf
end

function transform_stiffy(perm_Ks, wtf)
    num = length(perm_Ks)
    K_tran = zeros(12, 12, num)
    perm_K = zeros(12, 12, num)
    for i in 1:num
        perm_K[:,:,i] .= perm_Ks[i]
        K_tran[:,:,i] .= kron(I(4), wtf[:,:,i])
    end

    ## Using Einstein summation convention for passive transformation of stiffness matrix:
    #  D = Mᵗ(PK)M, where M = I₄ ⊗ T, and P is the permutation matrix acting on the original K
    @einsum D[j,k,i] := K_tran[l,j,i] * perm_K[l,m,i] * K_tran[m,k,i]
end

function build_big_stiffy(material, tubes, fem_mesh)
    num = length(tubes)

    # Building stiffness blocks
    Ks = [ tube_stiffness_matrix(material, [tube]) for tube in tubes ]
    # Ks[1]

    ## Testing simultaneous permutations of rows and columns:
    # 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
    # 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]
    inds    = [9,1,5,11,6,2,10,3,7,12,8,4]
    perm_Ks = [ K[inds,:][:,inds] for K in Ks ]

    # Compute transformation matrix
    wtf = axis_transformation(fem_mesh)

    # Transform the stiffness matrix
    D = transform_stiffy(perm_Ks, wtf)
end

## Stiffness matrix, loads, and constraints
D         = build_big_stiffy(aluminum, tubes, fem_mesh)
fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)

cons      = [span_num]
stiffness = build_stiffness_matrix(D, cons)

UnicodePlots.spy(stiffness)

## Solve system
displacements = solve_cantilever_beam(D, fem_loads, cons)

## Generate DataFrames
df_Fs = DataFrame(permutedims(fem_loads), :auto)
rename!(df_Fs, [:Fx, :Fy, :Fz, :Mx, :My, :Mz])

df_xs = DataFrame(permutedims(displacements), :auto)
rename!(df_xs, [:dx, :dy, :dz, :θx, :θy, :θz])

## Displacement transfer scheme
#==========================================================================================#

# I guess JJ doesn't know differential forms.
rotation_matrix(Ωx, Ωy, Ωz) = @SMatrix [  0  -Ωz  Ωy ;
                                          Ωz  0  -Ωx ;
                                         -Ωy  Ωx  0  ]

transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh) = permutedims(reduce(hcat, map(xyz -> xyz + dxs + Ts .* (xyz - fem_mesh), eachrow(vlm_mesh))))

##
θx, θy, θz = df_xs[!,4], df_xs[!,5], df_xs[!,6]
T          = rotation_matrix.(θx, θy, θz)
dxs        = SVector.(df_xs[!,1], df_xs[!,2], df_xs[!,3])
new_mesh   = transfer_displacements(dxs, T, vlm_mesh, fem_mesh)

## Plotting
#==========================================================================================#

using Plots
pyplot(dpi = 300)

# Beam circles
n_pts          = 20
circle3D(r, n) = [ (r*cos(θ), 0, r*sin(θ)) for θ in 0:2π/n:2π ];
circ           = circle3D(R, n_pts)

draw_tube(p1, p2, circ) = [ [ circ_pt .+ p1, circ_pt .+ p2 ] for circ_pt in circ ]

beam_pts   = zip(tupvector(fem_mesh[1:end-1]), tupvector(fem_mesh[2:end]))
left_pts   = [ draw_tube(pt[1], pt[2], circ) for pt in beam_pts ]

# Beam loads
fem_plot   = reduce(hcat, fem_mesh)
loads_plot = fem_loads

# Aerodynamic centers and forces
panel_plot = plot_panels(panels)
ac_plot    = reduce(hcat, vlm_acs)
force_plot = reduce(hcat, vlm_forces)

# Displacements
new_mesh_plot  = reduce(hcat, new_mesh)
new_panel_plot = plot_panels(make_panels(new_mesh)[:])

# Plot
b = aero_state.span_ref
plot(camera = (45, 45),
     xlim = (0, b/2),
    #  ylim = (-b/2, b/2),
     zlim = (0, b/2)
    )

# Panels
[ plot!(pans, color = :black, label = ifelse(i == 1, "Panels", :none)) for (i, pans) in enumerate(panel_plot) ]
[ plot!(pans, color = :brown, label = ifelse(i == 1, "Deflection", :none)) for (i, pans) in enumerate(new_panel_plot) ]

# Planform
plot!(wing_plan, color = :blue, label = "Planform")

# Beams
# [ plot!(reduce(vcat, pt), color = :green, label = ifelse(i == 1, "Beams", :none)) for (i, pt) in enumerate(left_pts) ]

# Forces
# quiver!(fem_plot[1,:], fem_plot[2,:], fem_plot[3,:],
#         quiver=(loads_plot[1,:], loads_plot[2,:], loads_plot[3,:] ) .* 0.1,
#         label = "Beam Forces")
# quiver!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:],
#         quiver=(force_plot[1,:], force_plot[2,:], force_plot[3,:]) .* 0.1,
#         label = "Panel Forces")
# scatter!(ac_plot[1,:], ac_plot[2,:], ac_plot[3,:], label = "Aerodynamic Centers")

# Beam axes
axes_plot  = reduce(hcat, (fem_mesh[1:end-1] + fem_mesh[2:end]) / 2)
ss         = @. normalize(fem_mesh[2:end] - fem_mesh[1:end-1])
ns         = ss .× fill(SVector(1,0,0), length(fem_mesh) - 1)
cs         = @. ss × ns
ss_plot    = reduce(hcat, ss)
cs_plot    = reduce(hcat, cs)
ns_plot    = reduce(hcat, ns)

quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:],
        quiver=(cs_plot[1,:], cs_plot[2,:], cs_plot[3,:]) .* 1e-1,
        color = :orange, label = :none)
quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:],
        quiver=(ss_plot[1,:], ss_plot[2,:], ss_plot[3,:]) .* 1e-1,
        color = :brown, label = :none)
quiver!(axes_plot[1,:], axes_plot[2,:], axes_plot[3,:],
        quiver=(ns_plot[1,:], ns_plot[2,:], ns_plot[3,:]) .* 1e-1,
        color = :red, label = :none)

# Displacements
scatter!(new_mesh_plot[1,:], new_mesh_plot[2,:], new_mesh_plot[3,:], label = "New Nodes")

plot!()