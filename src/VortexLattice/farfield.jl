## Farfield force evaluation
#==========================================================================================#

trefftz_potential(r_i :: SVector{3, <: Real}, r_j :: SVector{3, <: Real}, Γ_j :: Real, U_hat :: SVector{3, <: Real}) = let r = r_i .- r_j; Γ_j / 2π * U_hat × r / norm(r)^2 end

∇φ(r :: SVector{3, <: Real}, points, Γs :: AbstractVector{<: Real}, U_hat :: SVector{3, <: Real}) = sum(trefftz_potential.(Ref(r), points, Γs, Ref(U_hat)))

∂φ∂n(trefftz_line :: Line, points, Γs :: AbstractVector{<: Real}, normal, U_hat) = dot(∇φ(center(trefftz_line), points, Γs, U_hat), normal)

∂φ∂ns(trefftz_lines :: AbstractVector{Line}, Δφs :: AbstractVector{<: Real}, normals, U_hat) = ∂φ∂n.(trefftz_lines, (Ref ∘ points)(trefftz_lines), (Ref ∘ fwddiff)([ 0; -Δφs; 0 ]), normals, Ref(U_hat))

"""
    trefftz_forces(Γs, horseshoes, freestream, ρ)

Computes the aerodynamic forces in the Trefftz plane normal to the freestream given horseshoes, their associated strengths Γs, and a density ρ.
"""
function trefftz_forces(Γs, horseshoes :: AbstractArray{Horseshoe}, freestream :: Freestream, ρ :: Real)
    # Project trailing edge horseshoes' bound legs into Trefftz plane along wind axes   
    U_hat           = SVector(1, 0, 0)
    trefftz_lines   = body_to_wind_axes.(bound_leg.(horseshoes[end,:][:]), Ref(freestream))
    trefftz_vectors = vector.(trefftz_lines)
    trefftz_projs   = trefftz_vectors .- dot.(Ref(U_hat), trefftz_vectors) .* Ref(U_hat)
    normals         = Ref(U_hat) .× trefftz_projs

    # Compute directional derivatives of doublets in the normal direction
    @timeit "Sum Γs" Δφs = sum(Γs, dims = 1)[:]
    @timeit "∂φ/∂n" ∂φ_∂n = ∂φ∂ns(trefftz_lines, Δφs, normals, U_hat)

    # Compute forces    
    @timeit "Dihedrals" dihedrals = [ atan(vec[3], vec[2]) for vec in trefftz_projs ]
    @timeit "Projected Leg Norms" Δs = norm.(trefftz_projs)
    
    ΔφsΔs = Δφs .* Δs
    D_i   = -1/2 * ρ * sum(ΔφsΔs .* ∂φ_∂n)
    Y     = - ρ * freestream.mag * sum(ΔφsΔs .* sin.(dihedrals))
    L     = ρ * freestream.mag * sum(ΔφsΔs .* cos.(dihedrals))

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