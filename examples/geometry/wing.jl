## Wing example
using AeroFuse

wing = Wing(
    foils     = fill(naca4((2,4,1,2)), 3),  # Airfoils, type Foil
    chords    = [2.0, 1.6, 0.2],            # Chord lengths (m)
    twists    = [0.0, 0.0, 0.0],            # Twist angles (deg)
    spans     = [5., 0.6],                  # Span lengths (m)
    dihedrals = [5., 5.],                   # Dihedral angles (deg)
    sweeps    = [20.,20.],                  # Sweep angles (deg)
    w_sweep   = 0.25,                       # Normalized sweep location to chord ∈ [0,1]
    symmetry  = true,                       # Symmetry in x-z plane
    # flip      = true                      # Reflection about x-z plane
)

## Plotting
using Plots

plot(wing, aspect_ratio = 1, zlim = (-0.5, 0.5) .* span(wing), label = "Wing")