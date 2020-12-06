"""
Evaluates the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the global reference frame.
"""
stream_velocity(r, Ω :: SVector{3, Float64}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, V) = sum(velocity(r, horseshoe, Γ, V / norm(V)) for (horseshoe, Γ) ∈ zip(horseshoes, Γs)) .+ V .+ cross(Ω, r)

"""
Computes the streamlines from a given starting point, a Freestream, Horseshoes and their associated strengths Γs. The length of the streamline and the number of evaluation points must also be specified.
"""
function streamlines(point, freestream :: Freestream, Ω, horseshoes, Γs, length, num_steps)
    streamlines = fill(SVector{3, Float64}(0,0,0), num_steps)
    V = velocity(freestream)
    streamlines[1] = point
    cuck = x -> stream_velocity(x, Ω, horseshoes, Γs, V)
    @timeit "Iterating" for i ∈ 2:num_steps
        @timeit "Updating Velocity" update = cuck(streamlines[i-1])
        @timeit "Adding Streamline" streamlines[i] = streamlines[i-1] .+ (update / norm(update) * length / num_steps)
    end
    streamlines
end

"""
Computes the streamlines from the collocation points of given Horseshoes with the relevant previous inputs.
"""
streamlines(freestream :: Freestream, Ω, horseshoe_panels :: Array{<: Panel}, horseshoes :: Array{Horseshoe}, Γs :: Array{<: Real}, length :: Real, num_steps :: Integer) = [ streamlines(SVector(hs), freestream, Ω, horseshoes, Γs, length, num_steps) for hs ∈ collocation_point.(horseshoe_panels)[:] ]
