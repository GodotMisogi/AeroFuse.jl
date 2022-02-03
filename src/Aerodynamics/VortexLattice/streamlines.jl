"""
    stream_velocity(r, horseshoes, Γs, V, Ω)

Evaluate the total induced velocity at a point `r` given Horseshoes, vortex strengths `Γ`s, rotation rates `Ω`, and freestream flow vector `freestream` in the global reference frame.
"""
stream_velocity(r, horseshoes, Γs, V, Ω) = sum(x -> velocity(r, x[1], x[2], V / norm(V)), zip(horseshoes, Γs)) + V + Ω × r

"""
    streamlines(point, freestream :: Freestream, horseshoes, Γs, length, num_steps)

Compute the streamlines from a given starting point, a Freestream, Horseshoes and their associated strengths Γs with a specified length of the streamline and number of evaluation points.
"""
function streamlines(V, Ω, horseshoes, Γs, point, length, num_steps :: Integer)
    streamlines = fill(point, num_steps)
    f(x) = stream_velocity(x, horseshoes, Γs, V, Ω)
    for i ∈ 2:num_steps
        update = f(streamlines[i-1])
        streamlines[i] = streamlines[i-1] + (update / norm(update) * length / num_steps)
    end
    streamlines
end

streamlines(fs :: Freestream, refs :: References, vortices, Γs, points, length, num_steps :: Integer) = mapreduce(pt -> streamlines(refs.speed * velocity(fs), fs.omega, vortices, Γs, pt, length, num_steps), hcat, points)

streamlines(system :: VortexLatticeSystem, points, length, num_steps :: Integer) = streamlines(system.freestream, system.reference, system.vortices, system.circulations, points, length, num_steps)