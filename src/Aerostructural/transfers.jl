## Load transfer scheme
#==========================================================================================#

# Sum adjacent values
@views adjacent_adder(x1, x2) = [ [ x1[1] ]; x1[2:end] .+ x2[1:end-1]; [ x2[end] ] ]

# Compute moments for each section with local beam nodes as origins
section_moment(vlm_ac, fem_pts, half_vlm_force) = @. (vlm_ac - fem_pts) × half_vlm_force
@views section_moments(vlm_acs, fem_pts, half_vlm_forces) = sum(x -> section_moment(x[1], fem_pts, x[2]), zip(eachrow(vlm_acs), eachrow(half_vlm_forces)))

function compute_loads!(loads, vlm_acs, vlm_forces, fem_mesh)
    # Forces
    sec_forces = vec(sum(vlm_forces, dims = 1)) / 2
    beam_forces = adjacent_adder(sec_forces / 2, sec_forces / 2)

    # Moments
    M_ins = @views section_moments(vlm_acs, fem_mesh[1:end-1], vlm_forces / 2)
    M_outs = @views section_moments(vlm_acs, fem_mesh[2:end],   vlm_forces / 2)
    beam_moments = adjacent_adder(M_ins, M_outs)

    # Insert forces and moments into loads array
    loads[1:3,:] = combinedimsview(beam_forces)
    loads[4:6,:] = combinedimsview(beam_moments)

    nothing
end

function compute_loads(vlm_acs, vlm_forces, fem_mesh)
    T = promote_type(eltype(vlm_acs[1]), eltype(vlm_forces[1]), eltype(fem_mesh[1]))
    loads = MMatrix{6, length(fem_mesh), T}(undef)
    compute_loads!(loads, vlm_acs, vlm_forces, fem_mesh)

    return loads
end

# Generate load vector for FEM system
fem_load_vector(vlm_acs, vlm_forces, fem_mesh) = [ zeros(6); vec(compute_loads(vlm_acs, vlm_forces, fem_mesh)) ]


## Displacement transfer scheme
#==========================================================================================#

# Build cross product as an antisymmetric bilinear form. 
# (This is just a needlessly fancy way of saying antisymmetric matrix in 3 dimensions. The cross product, more generally, is actually the exterior product in higher dimensions.)
rotation_matrix(Ωx, Ωy, Ωz) = @SMatrix [  0  -Ωz  Ωy ;
                                          Ωz  0  -Ωx ;
                                         -Ωy  Ωx  0  ]

rotation_matrix(θs) = rotation_matrix.(θs[1,:], θs[2,:], θs[3,:])

# Transfer states by summing the displacements including rotations.
transfer_displacement(xyz, dx, rot, r) = xyz + dx + rot * (xyz - r)
transfer_displacements(dxs, Ts, chord_mesh, fem_mesh) = combinedimsview(map(xyz -> transfer_displacement.(xyz, dxs, Ts, fem_mesh), eachrow(chord_mesh)), (1))

@views mesh_translation(δs) = SVector.(δs[1,:], δs[2,:], δs[3,:])
@views mesh_rotation(δs) = rotation_matrix(δs[4:6,:])

# Make new horseshoes
function new_horseshoes(dxs, Ts, chord_mesh, camber_mesh, fem_mesh)
    @timeit "Chord Transfer" new_chord_mesh = transfer_displacements(dxs, Ts, chord_mesh, fem_mesh)
    @timeit "Camber Transfer" new_camber_mesh = transfer_displacements(dxs, Ts, camber_mesh, fem_mesh)
    @timeit "Chord Panels" cho_panels = make_panels(new_chord_mesh)
    @timeit "Normal Vectors" cam_panels = make_panels(new_camber_mesh)
    @timeit "Make Horseshoes" hs = map((cho, cam) -> Horseshoe(cho, normal_vector(cam)), cho_panels, cam_panels)

    return hs
end

new_horseshoes(Δs, chord_mesh, camber_mesh, fem_mesh) = new_horseshoes(mesh_translation(Δs), mesh_rotation(Δs), chord_mesh, camber_mesh, fem_mesh)