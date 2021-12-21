# Create the FEM mesh for the beam based on the VLM mesh
make_beam_mesh(vlm_mesh, fem_w) =  (1 - fem_w) * vlm_mesh[1,:] + fem_w * vlm_mesh[end,:]

# Axis transformation
function axis_transformation(stream, normie) 
    s_hat = normalize(stream)
    n_hat = normalize(normie)
    c_hat = s_hat × n_hat

    [ n_hat -c_hat -s_hat ]
end

# Compute local beam axes' transformation matrices using meshes (useless now?)
axis_transformation(fem_mesh :: Vector{SVector{3,T}}, vlm_mesh :: Array{SVector{3,T}}) where T <: Real = @. axis_transformation(fem_mesh[2:end] - fem_mesh[1:end-1], (vlm_mesh[end,2:end] - vlm_mesh[1,1:end-1]) × (vlm_mesh[1,2:end] - vlm_mesh[end,1:end-1]))

# Simultaneous permutations of rows and columns of the stiffness matrix
permute_stiffy(K) = let inds = [9,1,5,11,6,2,10,3,7,12,8,4]; @view K[inds,:][:,inds] end

# Index sets
# 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
# 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]

# Permute and (passively) transform stiffness matrix into global axes
transform_stiffy(K, axis) = let trans = kron(I(4), axis); permutedims(permutedims(trans) * permute_stiffy(K) * trans) end

# Simultaneously compute, permute and transform the stiffness matrices
build_big_stiffy(tubes, fem_mesh, vlm_mesh) = @views combinedimsview(@. transform_stiffy(tube_stiffness_matrix(tubes), axis_transformation(fem_mesh[2:end] - fem_mesh[1:end-1], (vlm_mesh[end,2:end] - vlm_mesh[1,1:end-1]) × (vlm_mesh[1,2:end] - vlm_mesh[end,1:end-1]))))

function evaluate_structures!(R_S, δ, xs, forces, fem_mesh, stiffness_matrix)
    # Build force vector with constraint for structures
    @timeit "FEM Loads" fem_loads = fem_load_vector(xs, forces, fem_mesh)

    # Structural residuals
    linear_residual!(R_S, stiffness_matrix, δ, fem_loads)

    fem_loads
end