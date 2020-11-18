## 
using Revise
includet("../src/FoilParametrization.jl")

foil = naca4((4,4,1,2))

num_secs = 3
foils = [ foil for i ∈ 1:num_secs ]
airfoils = Foil.(foils)

wing_chords = [0.18, 0.16, 0.08]
wing_twists = [2, 0, -2]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0, 11.3]
wing_sweeps = [1.14, 8]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
ρ = 1.225
ref = (0.25 * mean_aerodynamic_chord(wing), 0, 0)
uniform = Uniform(10.0, 5.0, -5.0)
@time horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 10, chord_num = 5, print = true);