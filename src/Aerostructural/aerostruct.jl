## Aerodynamic analysis
#==========================================================================================#

# lifting_line_forces(surf) = sum(surface_forces(surf), dims = 1)[:]
aerodynamic_forces(surfs) = reduce(vcat, (vec ∘ surface_forces).(surfs))

total_force(surfs) = sum(sum ∘ surface_forces, surfs) 

## Structural analysis
#==========================================================================================#

make_beam_mesh(vlm_mesh, fem_w = 0.35) =  (1 - fem_w) * vlm_mesh[1,:] + fem_w * vlm_mesh[end,:]

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

# Cantilever setup
function solve_beam_system(D, loads, constraint_indices)
    # Create the stiffness matrix from the array of individual stiffnesses
    # Also specifies the constraint location
    K = build_stiffness_matrix(D, constraint_indices)
    
    # Build force vector with constraint
    f = [ zeros(6); loads[:] ]

    # Solve FEM sys
    x = K \ f 

    # Throw away the junk values for the constraint
    reshape(x[7:end], 6, length(D[1,1,:]) + 1)
end

solve_beam_residual!(R, K, δ, F) = R .= K * δ - F

## Load-displacement transfer mechanisms
#==========================================================================================#

# Sum adjacent values
adjacent_joiner(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

section_moments(vlm_acs, fem_pts, half_vlm_forces) = sum(x -> (x[1] .- fem_pts) .× x[2], zip(eachrow(vlm_acs), eachrow(half_vlm_forces)))

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

# I guess JJ doesn't know differential forms.
rotation_matrix(Ωx, Ωy, Ωz) = @SMatrix [  0  -Ωz  Ωy ;
                                          Ωz  0  -Ωx ;
                                         -Ωy  Ωx  0  ]

rotation_matrix(θs) = rotation_matrix.(θs[1,:], θs[2,:], θs[3,:])

transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh) = permutedims(reduce(hcat, map(xyz -> xyz + dxs + Ts .* (xyz - fem_mesh), eachrow(vlm_mesh))))

transfer_normals(Ts, normals) = permutedims(reduce(hcat, map(normie -> (Ts[1:end-1] + Ts[2:end]) / 2 .* normie, eachrow(normals))))

## Weights and fuel loads
#==========================================================================================#

load_factor_residual(L, W, n) = L - n * W
load_factor_residual!(R, L, W, n) = R .= load_factor_residual(L, W, n)

## Coupled residual system
#==========================================================================================#

function solve_coupled_residual!(R, x, aero_system :: VLMSystem, aero_state :: VLMState, vlm_mesh, fem_mesh, weight, load_factor)
    n = (prod ∘ size)(horseshoes(aero_system))   # Get panel size

    # Unpack aerodynamic and structural variables
    Γ = @view x[1:n] 
    δ = @view x[n+1:end-1]

    # Get residual vector views
    R_A = @view R[1:n]
    R_S = @view R[n+1:end-1]
    R_W = @view R[end]
    
    # Compute displacements
    δs  = @views reshape(δ[7:end],  6, length(fem_mesh))
    dxs = @views SVector.(δs[1,:], δs[2,:], δs[3,:])
    Ts  = @views rotation_matrix.(δs[4,:], δs[5,:], δs[6,:])
    
    # Perturb VLM mesh and normals
    new_vlm_mesh   = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)

    # New VLM variables
    new_panels = make_panels(new_vlm_mesh)
    aero_system.surfaces[1].horseshoes = horseshoe_line.(new_panels)
    # aero_system.surfaces[1].normals    = transfer_normals(Ts, reshape(aero_system.surfaces[1].normals, size(new_panels)))

    # Aerodynamic residuals
    aero_state.alpha = x[end]
    solve_aerodynamic_residual!(R_A, Γ, aero_system, aero_surfs, aero_state)
    
    # Compute loads
    vlm_acs    = bound_leg_center.(horseshoes(aero_system.surfaces[1]))
    vlm_forces = surface_forces(aero_system.surfaces[1])

    # Compute FEM mesh
    fem_w     = 0.35
    fem_mesh  = make_beam_mesh(vlm_mesh, fem_w)
    Ls        = norm.(diff(fem_mesh)) # Beam lengths, m 
    tubes     = Tube.(Ref(aluminum), Ls, 0.1, 1e-3)


    # Create the stiffness matrix from the array of individual stiffnesses and specify constraints
    D    = build_big_stiffy(aluminum, tubes, fem_mesh)
    cons = [length(fem_mesh) ÷ 2]
    K    = build_stiffness_matrix(D, cons)
    
    # Build force vector with constraint
    fem_loads = compute_loads(vlm_acs, vlm_forces, fem_mesh)
    f = [ zeros(6); fem_loads[:] ]

    # Structural residuals
    solve_beam_residual!(R_S, K, δ, f)

    # Compute lift for load factor residual
    L = total_force(aero_surfs)[3]

    # Weight residual
    load_factor_residual!(R_W, L, weight, load_factor)

    R
end