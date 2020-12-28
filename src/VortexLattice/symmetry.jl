## Computations for symmetric cases
#==========================================================================================#

"""
    sym_influence_coefficient(r, horseshoe, panel_normal, V_hat, symmetry)

Computes the influence coefficient of the velocity of a Horseshoe with trailing lines in a given direction ``\\hat V`` at a point ``r`` projected to a normal vector, with symmetry in the ``x``-``z`` plane.
"""
function sym_influence_coefficient(r :: SVector{3, <: Real}, horseshoe :: Horseshoe, panel_normal :: SVector{3, <: Real}, V_hat :: SVector{3, <: Real})
    col_vel = velocity(r, horseshoe, 1., V_hat)
    mir_vel = (reflect_xz ∘ velocity)(reflect_xz(r), horseshoe, 1., V_hat)
    
    dot(col_vel .+ mir_vel, panel_normal)
end

"""
    sym_influence_matrix(colpoints, normals, horseshoes, V_hat)

Assembles the Aerodynamic Influence Coefficient (AIC) matrix given horseshoes, collocation points, associated normal vectors, a unit vector representing the freestream, with symmetry in the ``x``-``z`` plane.
"""
sym_influence_matrix(colpoints, normals, horseshoes :: AbstractVector{Horseshoe}, V_hat :: SVector{3, <: Real}) = @timeit "Matrix Construction" [ @timeit "Influence Coefficient" sym_influence_coefficient(r_i, horsie_j, n_i, V_hat) for (r_i, n_i) ∈ zip(colpoints, normals), horsie_j ∈ horseshoes ]


function sym_trailing_legs_velocities(r, line :: Line, Γ; direction = SVector(1., 0., 0.))
    a, b = r .- line.r1, r .- line.r2
    a_m, b_m = reflect_xz(a), reflect_xz(b)
   
    @timeit "Trailing Leg" trailing_legs_velocities(a, b, Γ, direction) .+ (reflect_xz ∘ trailing_legs_velocities)(a_m, b_m, Γ, direction)
end

function sym_trefftz_forces(Γs, horseshoes :: AbstractArray{Horseshoe}, freestream :: Freestream, ρ :: Real)
    # Symmetry condition
    horseshoes, Γs = symmetry ? [ horseshoes[:,end:1:-1] horseshoes ], [ Γs[:, end:1:-1] Γs ] : horseshoes, Γs

