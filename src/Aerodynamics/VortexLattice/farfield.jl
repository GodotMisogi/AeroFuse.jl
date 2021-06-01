## Farfield dynamics
#==========================================================================================#

trefftz_velocity(r_i, r_j, Γ_j) = let r = r_i - r_j; Γ_j / 2π * SVector(1, 0, 0) × r / norm(r)^2 end

trefftz_influence_matrix(centers, normals, points) = [ dot(trefftz_velocity(r_i, r_j2, 1.) - trefftz_velocity(r_i, r_j1, 1.), n_i) for (r_i, n_i) in zip(centers, normals), (r_j2, r_j1) in zip(points[2:end], points) ]

function doublet_normal_derivatives(trefftz_lines :: Vector{<: Line}, Δφs, normals)
	centers = center.(trefftz_lines)
	pts 	= points(trefftz_lines)
	AIC   	= trefftz_influence_matrix(centers, normals, pts) 
	∂φ_∂n 	= AIC * Δφs
end

function trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, α, β)
	U_ref     			= (Ref ∘ SVector)(1, 0, 0)

	# Transform to wind axes
	trefftz_lines       = @. body_to_wind_axes(bound_leg(horseshoes[end,:][:]), α, β)

	# Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes
	trefftz_vectors     = @. vector(trefftz_lines)
	trefftz_proj_vecs   = @. trefftz_vectors - dot(U_ref, trefftz_vectors) * U_ref

	# Normals, dihedral angles, and length
	normals             = @. normalize(U_ref × trefftz_proj_vecs)

	zs, ys 				= @. getindex(trefftz_proj_vecs, 3), getindex(trefftz_proj_vecs, 2)
	dihedrals 			= @. atan(zs, ys)
	Δs 					= norm.(trefftz_proj_vecs)

	trefftz_lines, normals, dihedrals, Δs
end

"""
	trefftz_compute(Δφs, Δs, ∂φ_∂n, θs, V, ρ, symmetry)

Compute the aerodynamic forces in the Trefftz plane given cumulative doublet strengths ``\\Delta \\phi``s, Trefftz panel lengths ``Δs``, doublet-normal directional derivatives ``\\partial \\phi / \\partial n``, Trefftz panel angles ``\\theta``s, the freestream speed and density ``V, \\rho``, and an option for symmetry.
"""
function trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, V, ρ)     
	D_i = - 1/2 * ρ * sum(@. Δφs * Δs * ∂φ_∂n)
	Y   = - ρ * V * sum(@. Δφs * Δs * sin(dihedrals))
	L   = ρ * V * sum(@. Δφs * Δs * cos(dihedrals))
	
	SVector(D_i, Y, L)
end

# ifelse(symmetry, SVector(D_i, 0, 2L), SVector(D_i, Y, L))

"""
	trefftz_forces(Γs, horseshoes, freestream, ρ)

Compute the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{<: Horseshoe}, speed, α, β, ρ)
	# Get projections of horseshoes into Trefftz plane with the associated normals, dihedral angles and lengths
	trefftz_lines, normals, dihedrals, Δs = trefftz_preprocessing(horseshoes, α, β)

	# Compute directional derivatives of doublets in the normal direction
	Δφs   = vec(sum(Γs, dims = 1))
	∂φ_∂n = doublet_normal_derivatives(trefftz_lines, Δφs, normals)

	# Compute forces
	trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, speed, ρ)
end