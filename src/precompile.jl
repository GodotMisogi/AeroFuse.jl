# Define wing
wing = Wing(
    foils     = [ naca4((0,0,1,2)) for i ∈ 1:2 ],
    chords    = [0.18, 0.16],
    twists    = [0., 0.],
    spans     = [0.25,],
    dihedrals = [5.],
    sweeps    = [1.14],
    symmetry  = true
)

# Get wing info
b = span(wing)
S = projected_area(wing)
c = mean_aerodynamic_chord(wing)
AR = aspect_ratio(wing)
λ = taper_ratio(wing)
mac = mean_aerodynamic_center(wing)


# Define freestream and reference values
fs = Freestream(2.0, 2.0, [0.0, 0.0, 0.0])
refs = References(
    speed    = 1.0, 
    area     = projected_area(wing), 
    span     = span(wing), 
    chord    = mean_aerodynamic_chord(wing), 
    density  = 1.225, 
    location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
)

aircraft = ComponentArray(wing = make_horseshoes(WingMesh(wing, [20], 5, span_spacing = [Sine(1); Sine()])))

# Evaluate stability case
system = VortexLatticeSystem(aircraft, fs, refs, false)
dv_data = freestream_derivatives(system)
dv_data = freestream_derivatives(system; axes = Wind())