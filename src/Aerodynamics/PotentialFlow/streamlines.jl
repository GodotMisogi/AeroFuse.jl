## Streamlines
#=================================================#

# Evaluate freestream, rotation, and solved-element induction in geometry axes.
stream_velocity(r, elements, strengths, V, Ω) = sum(x -> velocity(r, x[1], x[2], V / norm(V)), zip(elements, strengths)) + V + Ω × r

# Integrate the flow direction from a starting point with constant spatial steps.
function streamlines(V, Ω, elements, strengths, point, length, num_steps :: Integer)
    streamlines = fill(point, num_steps)
    f(x) = stream_velocity(x, elements, strengths, V, Ω)
    for i ∈ 2:num_steps
        update = f(streamlines[i-1])
        streamlines[i] = streamlines[i-1] + (update / norm(update) * length / num_steps)
    end
    streamlines
end

streamlines(fs :: Freestream, refs :: References, elements, strengths, points, length, num_steps :: Integer) = mapreduce(pt -> streamlines(refs.speed * velocity(fs), fs.omega, elements, strengths, pt, length, num_steps), hcat, points)

"""
    streamlines(system::PotentialFlowSystem, points, length, num_steps)

Integrate streamlines from `points` using freestream, rotation, and solved-element
induction in geometry axes. Prescribed fuselage thickness-source and propeller
slipstream fields are not included.
"""
streamlines(system :: PotentialFlowSystem, points, length, num_steps :: Integer) = streamlines(system.freestream, system.reference, system.elements, system.strengths, points, length, num_steps)
