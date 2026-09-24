# Scaling benchmark of the potential-flow solve and nearfield post-processing across backends:
# serial host (`backend = nothing`), multithreaded host (`CPU()`), and Apple GPU (`MetalBackend()`).
# Run from the repository root with threads enabled:
#   julia --project=test/gpu -t auto test/gpu/benchmark.jl

# %%
using AeroFuse
using ComponentArrays
using KernelAbstractions
using Metal
using Printf

wing = Wing(
    foils     = fill(naca4((0,0,1,2)), 2),
    chords    = [1.0, 0.6],
    twists    = [0.0, 0.0],
    spans     = [5.0] / 2,
    dihedrals = [11.39],
    sweeps    = [0.],
    symmetry  = true,
)

fs   = Freestream(alpha = 1.0, beta = 1.0, omega = zeros(3))
refs = References(
    speed    = 150.0,
    area     = projected_area(wing),
    span     = span(wing),
    chord    = mean_aerodynamic_chord(wing),
    density  = 1.225,
    location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.],
)

aircraft(n_span, n_chord) = ComponentVector(wing = elements(WingMesh(wing, [n_span], n_chord; span_spacing = Cosine()), Horseshoe()))

# Minimum wall time over `n` runs, in seconds.
besttime(f, n) = minimum(_ -> (t = time_ns(); f(); (time_ns() - t) / 1e9), 1:n)

function run_case(ac, backend; n = 3)
    sys = PotentialFlowSystem(ac, fs, refs; backend) # Warm-up (compilation) and reuse for post-processing
    t_solve = besttime(() -> PotentialFlowSystem(ac, fs, refs; backend), n)
    t_post  = besttime(() -> surface_coefficients(sys), n)
    return t_solve, t_post
end

# %%
backends = [("serial", nothing), ("CPU($(Threads.nthreads())t)", CPU()), ("Metal", MetalBackend())]
serial_limit = 6000 # Skip the serial path beyond this size to bound run time.

@printf("%8s  %-10s %10s %10s\n", "N", "backend", "solve [s]", "post [s]")
for (n_span, n_chord) in [(16, 10), (32, 20), (48, 30), (64, 40), (96, 50), (128, 60)]
    ac = aircraft(n_span, n_chord)
    N  = length(ac)
    for (label, backend) in backends
        isnothing(backend) && N > serial_limit && continue
        t_solve, t_post = run_case(ac, backend)
        @printf("%8d  %-10s %10.4f %10.4f\n", N, label, t_solve, t_post)
    end
end
