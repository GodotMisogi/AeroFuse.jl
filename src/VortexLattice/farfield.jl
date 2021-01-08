## Farfield dynamics
#==========================================================================================#

trefftz_potential(r_i :: SVector{3, <: Real}, r_j :: SVector{3, <: Real}, Γ_j :: Real) = let r = r_i - r_j; Γ_j / 2π * SVector(1, 0, 0) × r / norm(r)^2 end

∇φ(r :: SVector{3, <: Real}, points, Γs :: AbstractVector{<: Real}) = sum(trefftz_potential.(Ref(r), points, Γs))

∂φ∂n(trefftz_line :: Line, points, Γs :: AbstractVector{<: Real}, normal) = dot(∇φ(center(trefftz_line), points, Γs), normal)

∂φ∂ns(trefftz_lines :: AbstractVector{<: Line}, Δφs :: AbstractVector{<: Real}, normals) = ∂φ∂n.(trefftz_lines, (Ref ∘ points)(trefftz_lines), (Ref ∘ fwddiff)([ 0; -Δφs; 0 ]), normals)

function trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream)
    U_hat               = SVector(1, 0, 0)
    U_ref, free_ref     = Ref(U_hat), Ref(freestream)
    trefftz_lines       = @. body_to_wind_axes(bound_leg(horseshoes[end,:][:]), free_ref)

    trefftz_vectors     = @. vector(trefftz_lines)
    trefftz_proj_vecs   = @. trefftz_vectors - dot(U_ref, trefftz_vectors) * U_ref
    normals             = @. U_ref × trefftz_proj_vecs

    @timeit "Dihedrals" dihedrals = [ atan(vec[3], vec[2]) for vec in trefftz_proj_vecs ]
    @timeit "Projected Leg Norms" Δs = norm.(trefftz_proj_vecs)

    trefftz_lines, trefftz_proj_vecs, normals, dihedrals, Δs
end

function trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, V, ρ, symmetry)     
    D_i = - 1/2 * ρ * sum(@. Δφs * Δs * ∂φ_∂n)
    Y   = - ρ * V * sum(@. Δφs * Δs * sin(dihedrals))
    L   = ρ * V * sum(@. Δφs * Δs * cos(dihedrals))
    
    symmetry ? SVector(D_i, 0, 2L) : SVector(D_i, Y, L)
end

"""
    trefftz_forces(Γs, horseshoes, freestream, ρ)

Computes the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream, ρ :: Real, symmetry :: Bool)
    # Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes
    trefftz_lines, trefftz_proj_vecs, normals, dihedrals, Δs = trefftz_preprocessing(horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream)

    # Compute directional derivatives of doublets in the normal direction
    @timeit "Sum Γs" Δφs = vec(sum(Γs, dims = 1))
    @timeit "∂φ/∂n" ∂φ_∂n = ∂φ∂ns(trefftz_lines, Δφs, normals)

    # Compute forces
    trefftz_compute(Δφs, Δs, ∂φ_∂n, dihedrals, freestream.mag, ρ, symmetry)
end 

"""
    trefftz_forces(Γs, horseshoes, freestream, r_ref, ρ)

Compute farfield forces and moments in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, a reference point for calculation of moments, and a density ρ.
"""
function farfield_dynamics(Γs :: AbstractArray{<: Real}, horseshoes :: AbstractArray{<: Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225, symmetry = false)
    @timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, freestream, ρ, symmetry)
    @timeit "Trefftz Moment" trefftz_moment = r_ref × trefftz_force

    trefftz_force, trefftz_moment
end