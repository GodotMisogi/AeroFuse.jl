#---------------------------Farfield evaluations---------------------------#

trefftz_potential(r_i :: SVector{3, <: Real}, r_j :: SVector{3, <: Real}, Γ_j :: Real, U_hat :: SVector{3, <: Real}) = let r = r_i .- r_j; Γ_j / 2π * U_hat × r / dot(r, r) end

∇φ(r :: SVector{3, <: Real}, Γ :: Real, trefftz_lines, U_hat :: SVector{3, <: Real}) = sum(trefftz_potential(r, r_jp½, Γ, U_hat) for r_jp½ ∈ point2.(trefftz_lines))

∂φ∂n(trefftz_lines, Γs, normals, U_hat) = dot.((∇φ(center(rline_i), Γ_i, trefftz_lines, U_hat) for (rline_i, Γ_i) ∈ zip(trefftz_lines, Γs)), normals)

"""
    trefftz_forces(Γs, horseshoes, freestream, ρ)

Computes the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{Horseshoe}, freestream :: Freestream, ρ :: Real)

    # Project horseshoes' bound legs into Trefftz plane along freestream
    @timeit "Sum Γs" Δφs = sum(Γs, dims = 1)
    trefftz_lines = bound_leg.(horseshoes[end,:])[:]
    trefftz_vectors = vector.(trefftz_lines)

    U = velocity(freestream)
    Us = fill(U, length(trefftz_lines))
    normals = normalize.(Us .× trefftz_vectors)

    # Compute directional derivatives of doublets
    @timeit "∂φ/∂n" ∂φ_∂n = ∂φ∂n(trefftz_lines, Δφs, normals, normalize(U))

    # Compute forces    
    @timeit "Dihedrals" dihedrals = [ atan(vec[3], vec[2]) for vec in trefftz_vectors ]
    @timeit "Projected Leg Norms" Δs = norm.(trefftz_vectors)
    
    pots_lens = Δφs .* Δs
    D_i = -0.5 * ρ * sum(pots_lens .* ∂φ_∂n)
    Y = - ρ * sum(pots_lens .* sin.(dihedrals))
    L = ρ * sum(pots_lens .* cos.(dihedrals))

    SVector(D_i, Y, L)
end 

"""
    trefftz_forces(Γs, horseshoes, freestream, r_ref, ρ)

Compute farfield forces and moments in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, a reference point for calculation of moments, and a density ρ.
"""
function farfield_dynamics(Γs :: AbstractArray{<: Real}, horseshoes :: AbstractArray{Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225)
    @timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, freestream, ρ)
    @timeit "Trefftz Moment" trefftz_moment = r_ref × trefftz_force

    trefftz_force, trefftz_moment
end