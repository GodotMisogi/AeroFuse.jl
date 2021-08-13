# Create the FEM mesh for the beam based on the VLM mesh
make_beam_mesh(vlm_mesh, fem_w) =  (1 - fem_w) * vlm_mesh[1,:] + fem_w * vlm_mesh[end,:]

# Axis transformation of stiffness matrices (NEEDS IMPROVEMENT)
function axis_transformation(fem_mesh, vlm_mesh)
    # Compute local beam axes
    ss = @. normalize(fem_mesh[2:end] - fem_mesh[1:end-1])
    ns = @. normalize((vlm_mesh[end,2:end] - vlm_mesh[1,1:end-1]) × (vlm_mesh[1,2:end] - vlm_mesh[end,1:end-1]))
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

# Using Einstein summation convention for passive transformation of stiffness matrix: (NEEDS IMPROVEMENT) 
# D = Mᵗ(PK)M, where M = I₄ ⊗ T, and P is the permutation matrix acting on the original K
function transform_stiffy(perm_Ks, wtf)
    num = length(perm_Ks)
    K_tran = zeros(12, 12, num)
    perm_K = zeros(12, 12, num)
    for i in 1:num
        perm_K[:,:,i] .= perm_Ks[i]
        K_tran[:,:,i] .= kron(I(4), wtf[:,:,i])
    end
    @einsum D[j,k,i] := K_tran[l,j,i] * perm_K[l,m,i] * K_tran[m,k,i]
end

# Simultaneous permutations of rows and columns of the stiffness matrix:
# 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
# 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]
function permute_stiffy(Ks) 
    inds    = [9,1,5,11,6,2,10,3,7,12,8,4]
    perm_Ks = [ K[inds,:][:,inds] for K in Ks ] 
end

# Build the full stiffness matrix with permutations and axis transformations
function build_big_stiffy(tubes, fem_mesh, vlm_mesh)
    # Building stiffness blocks
    Ks = [ tube_stiffness_matrix(tube) for tube in tubes ]

    perm_Ks = permute_stiffy(Ks)

    # Compute transformation matrix
    wtf = axis_transformation(fem_mesh, vlm_mesh)

    # Transform the stiffness matrix
    transform_stiffy(perm_Ks, wtf)    
end