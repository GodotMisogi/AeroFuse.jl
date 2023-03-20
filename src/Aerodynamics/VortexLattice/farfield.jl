## Farfield dynamics
#==========================================================================================#

"""
    farfield_velocity(r_i, r_j, Γ_j)

Compute the induced velocity at point ``r_i`` of the wake in the Trefftz plane due to the strength ``Γ_j`` at ``r_j``. 
"""
function farfield_velocity(r_i, r_j, Γ_j) 
    T = promote_type(eltype(r_i), eltype(r_j), typeof(Γ_j))
    r = r_i - r_j
    
    Γ_j / 2π * SVector{3,T}(1, 0, 0) × r / norm(r)^2 
end

"""
    farfield_influence_matrix(centers, normals, points)

Compute the aerodynamic influence coefficient matrix of the wake in the Trefftz plane given the center points, normal vectors, and the points of the wake panels.
"""
@views farfield_influence_matrix(centers, normals, points) = [ dot(farfield_velocity(r_i, r_j2, 1.) - farfield_velocity(r_i, r_j1, 1.), n_i) for (r_i, n_i) in zip(centers, normals), (r_j2, r_j1) in zip(points[2:end], points) ]

"""
    doublet_normal_derivatives(wake_points, Δφs, normals)

Compute the normal derivative strengths of the doublets given the wake points, net doublet strengths ``Δφ``s, and associated normal vectors.
"""
@views function doublet_normal_derivatives(wake_points, Δφs, normals)
    centers = @. (wake_points[1:end-1] + wake_points[2:end]) / 2
    AIC     = farfield_influence_matrix(centers, normals, wake_points)
    ∂φ_∂n   = AIC * Δφs

    return ∂φ_∂n
end

project_vector(vector, U) = vector - dot(U, vector) * U
normal(wake_proj_vector, U) = normalize(U × wake_proj_vector)
@views dihedral(wake_proj_vec) = atan(wake_proj_vec[3], wake_proj_vec[2])

"""
    trefftz_plane_quantities(vortices, α, β)

Project `Horseshoe`s into the Trefftz plane aligned with the wind axes angles ``α,~β``, and compute normal vectors, projection angles and lengths.
"""
function trefftz_plane_quantities(vortices, α, β)
    # Reference velocity for broadcasting
    U_ref = (Ref ∘ SVector)(1, 0, 0)

    # Transform wake to wind axes
    wake_horsies = @views vec(vortices[end,:])                  # Trailing edge
    wake_points  = [ r1.(wake_horsies); [r2(wake_horsies[end])] ] # Trailing edge points
    wake_wind    = geometry_to_wind_axes.(wake_points, α, β)      # Wind axis transformation

    # Project trailing edge vortices' bound legs into Trefftz plane along wind axes
    wake_vectors   = @. wake_wind[2:end] - wake_wind[1:end-1] # Vectors
    wake_proj_vecs = @. project_vector(wake_vectors, U_ref)   # Project to wind axes

    # Normals, dihedral angles, and projection lengths
    wake_normals   = @. normal(wake_proj_vecs, U_ref)
    wake_dihedrals = @. dihedral(wake_proj_vecs)
    wake_lengths   = @. norm(wake_proj_vecs)

    return wake_points, wake_normals, wake_dihedrals, wake_lengths
end

"""
    compute_farfield_forces(Δφs, Δs, ∂φ_∂n, θs, V, ρ)

Compute the aerodynamic forces in the Trefftz plane given cumulative doublet strengths ``Δφ``s, Trefftz panel lengths ``Δs``, doublet-normal directional derivatives ``∂φ/∂n``, Trefftz panel angles ``θ``s, the freestream speed ``V`` and density ``ρ``.
"""
function compute_farfield_forces(Δφs, Δs, ∂φ_∂n, θs, V, ρ)
    D_i = - 1/2 * ρ * sum(@. Δφs * Δs * ∂φ_∂n)
    Y = - ρ * V * sum(@. Δφs * Δs * sin(θs))
    L = ρ * V * sum(@. Δφs * Δs * cos(θs))
    
    SVector(D_i, Y, L)
end

"""
    farfield_forces(Γs, horseshoes, freestream, ρ)

Compute the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths ``Γ``s, and a density ``ρ``.
"""
function farfield_forces(Γs, horseshoes, speed, α, β, ρ)
    # Get projections of horseshoes into Trefftz plane with the associated normals, dihedral angles and lengths
    wake_points, normals, dihedrals, Δs = trefftz_plane_quantities(horseshoes, α, β)

    # Compute directional derivatives of doublets in the normal direction
    Δφs = vec(sum(Γs, dims = 1))
    ∂φ_∂n = doublet_normal_derivatives(wake_points, Δφs, normals)

    # Compute forces
    return compute_farfield_forces(Δφs, Δs, ∂φ_∂n, dihedrals, speed, ρ)
end