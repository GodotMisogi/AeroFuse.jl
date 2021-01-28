## 
using Revise
using StaticArrays
using AeroMDAO

## Wing section setup
foil        =   naca4((2,4,1,2))
wing_right  =   HalfWing(Foil.(foil for i ∈ 1:3),   # Foils
                        [0.18, 0.16, 0.08],         # Chords
                        [2., 0., -2.],              # Twists
                        [0.5, 0.5],                 # Spans
                        [0., 11.3],                 # Dihedrals
                        [1.14, 8.])                 # Sweeps

wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
ρ = 1.225
ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
Ω = [0.0, 0.0, 0.0]
uniform = Freestream(10.0, 5.0, 0.0, Ω)
nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ρ, ref; span_num = 10, chord_num = 10) 

begin
    println("\nNearfield:")
    print_dynamics(nf_coeffs...)
    println("\nFarfield:")
    print_dynamics(ff_coeffs...)
end

##
using PlotlyJS

# num_points = 50
# max_z = 0.1
# y = span(wing) / 2 - 0.05
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

span_points = 20
init        = trailing_chopper(wing.right, span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz])  ; 
                init .+ Ref([dx, dy, -dz]) ];

##
layout, plt = plot_case(horseshoe_panels, camber_panels, Γs, horseshoes, uniform, seed, 5., 100);
plt
