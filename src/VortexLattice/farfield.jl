
#---------------------------Farfield evaluations---------------------------#

trefftz_potential(r_i, r_j, Γ_j, Û) = let r = r_i .- r_j; Γ_j/2π * cross(Û, r) / dot(r, r) end

trefftz_matrix(trefftz_lines, normals, Û) = [ dot(trefftz_potential(center(tline_i), tline_j.r1, 1., Û), n̂_i) for (tline_i, n̂_i) in zip(trefftz_lines, normals), tline_j in trefftz_lines ]

project_yz(line :: Line) = Line(SVector(0, line.r1[2], line.r1[3]), SVector(0, line.r2[2], line.r2[3]))

body_to_wind_axes(line :: Line, freestream :: Freestream) = Line(body_to_wind_axes(line.r1, freestream), body_to_wind_axes(line.r2, freestream)) 

"""
Computes the aerodynamic forces in the Trefftz plane normal to the freestream.
"""
function trefftz_forces(Γs, horseshoes :: Array{Horseshoe}, freestream :: Freestream, ρ)

    # lines = bound_leg_vector.(horseshoes[end,:])

    # projs = dot.(vels, lines)
    # projected_vecs = lines .- (projs ./ norm.(projs) .* vels

    # lines = [ body_to_wind_axes(horseshoe.bound_leg, freestream) for horseshoe in horseshoes[end,:] ][:]
    # trefftz_lines = project_yz.()

    # Project horseshoes' bound legs into Trefftz plane along freestream
    U = velocity(freestream)
    trefftz_lines = [ horseshoe.bound_leg for horseshoe in horseshoes[end,:] ][:]
    trefftz_vectors = vector.(trefftz_lines)

    Us = repeat([U], length(trefftz_lines))
    normals = cross.(Us, trefftz_vectors)
    normals = normals ./ norm.(normals)

    # Compute matrices
    @timeit "Dihedrals" dihedrals = [ atan(vec[3], vec[2]) for vec in trefftz_vectors ]
    @timeit "Projected Leg Norms" Δs = norm.(trefftz_vectors)
    @timeit "Sum Γs" Δφs = sum(Γs, dims = 1)
    @timeit "Trefftz AIC" AIC = trefftz_matrix(trefftz_lines, normals, normalize(U))

    # Solve system
    @timeit "∂φ/∂n" ∂φ_∂n = AIC * Δφs[:]

    # Compute forces
    pots_lens = Δφs .* Δs
    D_i = -0.5 * ρ * sum(pots_lens .* ∂φ_∂n)
    Y = - ρ * sum(pots_lens .* sin.(dihedrals))
    L = ρ * sum(pots_lens .* cos.(dihedrals))

    SVector(D_i, Y, L)
end 

"""
Compute farfield forces and moments.
"""
function farfield_dynamics(Γs :: Array{<: Real}, horseshoes :: Array{Horseshoe}, freestream :: Freestream, r_ref, ρ = 1.225)
    @timeit "Trefftz Force" trefftz_force = trefftz_forces(Γs, horseshoes, freestream, ρ)
    @timeit "Trefftz Moment" trefftz_moment = cross(r_ref, trefftz_force)

    trefftz_force, trefftz_moment
end