## Aircraft geometry
using AeroFuse

## Surfaces

# Wing
wing = Wing(
    foils     = fill(naca4(2,4,1,2), 2),
    chords    = [1.0, 0.6],
    twists    = [2.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.],
    sweeps    = [5.],
    symmetry  = true,
    # flip      = true
);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25],
    dihedrals = [0.],
    sweeps    = [6.39],
    position  = [4., 0, -0.1],
    angle     = -2.,
    axis      = [0., 1., 0.],
    symmetry  = true
)

# Vertical tail
vtail = Wing(
    foils     = fill(naca4(0,0,0,9), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.0],
    dihedrals = [0.],
    sweeps    = [7.97],
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.]
)

# Fuselage
l_fuselage = 8.      # Length (m)
h_fuselage = 0.6     # Height (m)
w_fuselage = 0.7     # Width (m)

# Chordwise locations and corresponding radii
lens = [0.0, 0.005, 0.01, 0.03, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.98, 1.0]
rads = [0.05, 0.15, 0.25, 0.4, 0.8, 1., 1., 1., 1., 0.85, 0.3, 0.01] * w_fuselage / 2

fuse = Fuselage(l_fuselage, lens, rads, [-3.0, 0., 0.])


## Plotting
using Plots
plotlyjs()

plot(
    aspect_ratio = 1,
    zlim = (-0.5, 0.5) .* span(wing)
)
plot!(wing, label = "Wing")
plot!(htail, label = "Horizontal Tail")
plot!(vtail, label = "Vertical Tail")
plot!(fuse, color = :grey, label = "Fuselage")