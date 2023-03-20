# Wing
wing = Wing(
    foils     = fill(naca4((0,0,1,2)), 2),
    chords    = [1.0, 0.6],
    twists    = [0.0, 0.0],
    spans     = [5.0] / 2,
    dihedrals = [11.39],
    sweeps    = [0.],
    symmetry  = true, 
);

# Horizontal tail
htail = Wing(
    foils     = fill(naca4((0,0,1,2)), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25] / 2,
    dihedrals = [0.],
    sweeps    = [6.39],
    position  = [4., 0, 0],
    angle     = -2.,
    axis      = [0., 1., 0.],
    symmetry  = true
)

# Vertical tail
vtail = Wing(
    foils     = fill(naca4((0,0,0,9)), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.0],
    dihedrals = [0.],
    sweeps    = [7.97],
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.],
)

# Assembly
wing_mesh = WingMesh(wing, [32], 10; span_spacing = Cosine())
htail_mesh = WingMesh(htail, [12], 6; span_spacing = Cosine())
vtail_mesh = WingMesh(vtail, [5], 6; span_spacing = Cosine())

## Reference quantities
fs = Freestream(
    alpha    = 1.0, 
    beta     = 1.0, 
    omega    = zeros(3)
)

refs = References(
    speed    = 150.0,
    area     = projected_area(wing),
    span     = span(wing),
    chord    = mean_aerodynamic_chord(wing),
    density  = 1.225,
    location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
)
