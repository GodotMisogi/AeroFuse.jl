## Wing example
using AeroFuse

wing = Wing(
    foils = fill(naca4((2, 4, 1, 2)), 3),   # Airfoils, type Foil
    chords = [2.0, 1.6, 0.8],                # Chord lengths (m)
    twists = [0.0, 0.0, 0.0],                # Twist angles (deg)
    spans = [3.0, 1.2],                     # Span lengths (m)
    dihedrals = [5.0, 5.0],                     # Dihedral angles (deg)
    sweeps = [20.0, 40.0],                   # Sweep angles (deg)
    # controls  = [Flap(0, 0.75), Aileron(-30, 0.75)],
    sweep_ratio = 0.25,                         # Chord length fraction of sweep location
    symmetry = true,                            # Symmetry in x-z plane
    # flip      = true                          # Reflection about x-z plane
)

## Create symmetric wing instead
using Accessors

wing = @set wing.symmetry = true

## Evaluate geometric properties
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing) # Aerodynamic center
S_w = projected_area(wing) # Projected area
b_w = span(wing) # Wingspan
c_w = mean_aerodynamic_chord(wing) # Mean aerodynamic chord
tau_w = taper_ratio(wing) # Taper ratio

lambda_w_c4 = sweeps(wing, 0.25) # Quarter-chord sweep angles
lambda_w_c2 = sweeps(wing, 0.5) # Half-chord sweep angles

ct = camber_thickness(wing, 60) # Camber-thickness distribution
coords = coordinates(wing) # Leading and trailing edge coordinates

## Plotting
using Plots
plt = plot(
    size = (800, 600),
    aspect_ratio = 1,
    zlim = (-0.5, 0.5) .* span(wing),
    camera = (30, 60),
)
plot!(
    wing,
    label = "Wing",
    # mac = false # Disable mean aerodynamic center plot
)

# savefig(plt, "plots/wing_geom.pdf")

##