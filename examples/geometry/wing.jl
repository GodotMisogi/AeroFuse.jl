## Wing example
using AeroFuse

wing = Wing(
    foils     = fill(naca4((2,4,1,2)), 3),
    chords    = [2.0, 1.6, 0.2],
    twists    = [0.0, 0.0, 0.0],
    spans     = [5., 0.6],
    dihedrals = [5., 5.],
    sweeps    = [20.,20.],
    w_sweep   = 0.25, # Quarter-chord sweep
    symmetry  = true,
    # flip      = true
)

## Plotting
using Plots

plot(wing, aspect_ratio = 1, zlim = (-0.5, 0.5) .* span(wing), label = "Wing")