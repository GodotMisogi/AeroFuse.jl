## Wing analysis case
using AeroFuse
using Plots
using StructArrays

## Wing section setup
wing = Wing(
    foils = fill(naca4((0, 0, 1, 2)), 2),
    chords = [2.2, 1.8],
    twists = [0.0, 0.0],
    spans = [7.5],
    dihedrals = [0.0],
    sweeps = [3.0528],
    sweep_ratio = 0.0, # Quarter-chord sweep
    symmetry = true,
)

wing_mesh = WingMesh(wing, [24], 8, span_spacing = Uniform());

# Freestream conditions
fs = Freestream(
    alpha = 2.0, # deg
    beta = 2.0, # deg
    omega = [0.0, 0.0, 0.0],
)

# Reference values
ref = References(
    speed = 1.0, # m/s
    density = 1.225, # kg/m³
    viscosity = 1.5e-5, # ???
    area = projected_area(wing), # m²
    span = span(wing), # m
    chord = mean_aerodynamic_chord(wing), # m
    location = mean_aerodynamic_center(wing), # m
)

## Horseshoes
ac_hs = ComponentVector(wing = make_horseshoes(wing_mesh))
system = VortexLatticeSystem(ac_hs, fs, ref, false)
print_coefficients(nearfield(system), farfield(system))

## Vortex rings
ac_vs = ComponentVector(wing = make_vortex_rings(wing_mesh))
sys = VortexLatticeSystem(ac_vs, fs, ref, false)
print_coefficients(nearfield(sys), farfield(sys))
