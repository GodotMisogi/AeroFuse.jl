## 
# using BenchmarkTools
# using ProfileView
using Revise
using StaticArrays
using TimerOutputs
using AeroMDAO

## Wing section setup
foil = naca4((4,4,1,2))
wing_right  = HalfWing(Foil.(foil for i ∈ 1:3), # Foils
                        [0.18, 0.16, 0.08],     # Chords
                        [2., 0., -2.],          # Twists
                        [0.5, 0.5],             # Spans
                        [0., 11.3],             # Dihedrals
                        [1.14, 8.])             # Sweeps
wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
# reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 5.0, 0.0, Ω)
@time coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10, print = true) 

# print_timer();

##
using PlotlyJS

##
layout = Layout(title = "Vortex Lattice",
                scene=attr(aspectratio=attr(x=1,y=1,z=1)),
                )

## Streamlines
reset_timer!()

@timeit "Computing Streamlines" trace_streams = trace_streamlines(uniform, horseshoe_panels[:], horseshoes[:], Γs[:], 2, 100);

print_timer()

##
trace_horses = trace_panels(horseshoe_panels[:], Γs[:])
trace_cambers = trace_panels(camber_panels[:])
trace_wing = trace_surface(wing)

PlotlyJS.plot([ 
                  [ trace for trace in trace_horses ]...,
                  [ trace for trace in trace_cambers ]...,
                  [ trace for trace in trace_streams ]...,
                  # [ trace for trace in trace_wing ]...,
              ],
              layout)