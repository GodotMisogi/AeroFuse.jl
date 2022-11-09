# Create the FEM mesh for the beam based on the VLM mesh
@views make_beam_mesh(chord_mesh, fem_w) = (1 - fem_w) * chord_mesh[1,:] + fem_w * chord_mesh[end,:]

# Axis transformation
function axis_transformation!(mat, stream, normie) 
    s_hat = normalize(stream)
    n_hat = normalize(normie)
    c_hat = s_hat × n_hat

    mat[:,1] = n_hat
    mat[:,2] = -c_hat
    mat[:,3] = -s_hat

    nothing
end

function axis_transformation(stream, normie) 
    T = promote_type(eltype(stream), eltype(normie))
    mat = zeros(T, 3, 3)

    axis_transformation!(mat, stream, normie)

    return mat
end


# Compute local beam axes' transformation matrices using meshes (useless now?)
@views axis_transformation(fem_mesh :: Vector{SVector{3,T}}, chord_mesh :: Array{SVector{3,T}}) where T <: Real = @. axis_transformation(fem_mesh[2:end] - fem_mesh[1:end-1], (chord_mesh[end,2:end] - chord_mesh[1,1:end-1]) × (chord_mesh[1,2:end] - chord_mesh[end,1:end-1]))

# Simultaneous permutations of rows and columns of the stiffness matrix
permute_stiffy(K) = let inds = [9,1,5,11,6,2,10,3,7,12,8,4]; @view K[inds,:][:,inds] end

# Index sets
# 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
# 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]

# Permute and (passively) transform stiffness matrix into global axes
transform_stiffy(K, axis) = let trans = kron(I(4), axis); permutedims(permutedims(trans) * permute_stiffy(K) * trans) end

# Simultaneously compute, permute and transform the stiffness matrices
@views build_big_stiffy(tubes, fem_mesh, chord_mesh) = combinedimsview(@. transform_stiffy(tube_stiffness_matrix(tubes), axis_transformation(fem_mesh[2:end] - fem_mesh[1:end-1], (chord_mesh[end,2:end] - chord_mesh[1,1:end-1]) × (chord_mesh[1,2:end] - chord_mesh[end,1:end-1]))))