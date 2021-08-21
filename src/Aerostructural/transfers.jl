## Load transfer scheme
#==========================================================================================#

# Sum adjacent values
adjacent_adder(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

# Compute moments for each section with local beam nodes as origins
section_moments(vlm_acs, fem_pts, half_vlm_forces) = sum(x -> (x[1] .- fem_pts) .× x[2], zip(eachrow(vlm_acs), eachrow(half_vlm_forces)))

function compute_loads(vlm_acs, vlm_forces, fem_mesh)
    # Forces
    sec_forces   = sum(vlm_forces, dims = 1)[:] / 2
    beam_forces  = adjacent_adder(sec_forces / 2, sec_forces / 2)

    # Moments
    M_ins        = @views section_moments(vlm_acs, fem_mesh[1:end-1], vlm_forces / 2)
    M_outs       = @views section_moments(vlm_acs, fem_mesh[2:end],   vlm_forces / 2)
    beam_moments = adjacent_adder(M_ins, M_outs)

    # Concatenate forces and moments into loads array
    [ reduce(hcat, beam_forces); reduce(hcat, beam_moments) ]
end

# Generate load vector for FEM system
fem_load_vector(vlm_acs, vlm_forces, fem_mesh) = [ zeros(6); compute_loads(vlm_acs, vlm_forces, fem_mesh)[:] ]


## Displacement transfer scheme
#==========================================================================================#

# Build cross product as an antisymmetric bilinear form. 
# (This is just a fancy way of saying antisymmetric matrix in 3 dimensions. The cross product, more generally, is actually the exterior product in higher dimensions.)
rotation_matrix(Ωx, Ωy, Ωz) = @SMatrix [  0  -Ωz  Ωy ;
                                          Ωz  0  -Ωx ;
                                         -Ωy  Ωx  0  ]

rotation_matrix(θs) = rotation_matrix.(θs[1,:], θs[2,:], θs[3,:])

# Transfer states by summing the displacements and the cross product of the rotations with the .
transfer_displacement(xyz, dx, rot, r) = xyz + dx + rot * (xyz - r)
transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh) = permutedims(reduce(hcat, map(xyz -> transfer_displacement.(xyz, dxs, Ts, fem_mesh), eachrow(vlm_mesh))))

translations_and_rotations(δs) = @views SVector.(δs[1,:], δs[2,:], δs[3,:]), rotation_matrix(δs[4:6,:])

# Make new horseshoes
function new_horseshoes(dxs, Ts, vlm_mesh, cam_mesh, fem_mesh)
    new_vlm_mesh = transfer_displacements(dxs, Ts, vlm_mesh, fem_mesh)
    new_cam_mesh = transfer_displacements(dxs, Ts, cam_mesh, fem_mesh)
    Horseshoe.(make_panels(new_vlm_mesh), panel_normal.(make_panels(new_cam_mesh)))
end