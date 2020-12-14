## 
using Revise
using StaticArrays
using BenchmarkTools
using TimerOutputs
using ProfileView
using AeroMDAO

## Wing section setup
foil = naca4((4,4,1,2))
num_secs = 3
foils = [ foil for i ∈ 1:num_secs ]

airfoils = Foil.(foils)
wing_chords = [0.18, 0.16, 0.08]
wing_twists = [2., 0., -2.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0., 11.3]
wing_sweeps = [1.14, 8.]

# wing_right = HalfWing(airfoils, wing_chords, wing_twists, wing_spans, wing_dihedrals, wing_sweeps)
wing_right = HalfWing(Foil.(naca4((2,4,1,2)) for i ∈ 1:5),
                      [0.0639628599561049, 0.06200381820887121, 0.05653644812231768, 0.04311297779068357, 0.031463501535620116],
                      [0., 0., 0., 0., 0.],
                      [0.2, 0.2, 0.2, 0.2],
                      [0., 0., 0., 0.],
                      [0., 0., 0., 0.])
wing = Wing(wing_right, wing_right)
print_info(wing)

## Assembly
reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 5.0, 0.0, Ω)
@time coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10, print = true) 

print_timer();

##
using PlotlyJS

##
layout = Layout(
                title = "Vortex Lattice",
                scene=attr(aspectmode="manual", aspectratio=attr(x=1,y=1,z=1)),
                zlim=(-0.1, 5.0)
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
        #       [ trace for trace in trace_wing ]...,
              ],
              layout)

##
# using Plots
# plotlyjs()

# ##
# plot(xaxis = "x", yaxis = "y", zaxis = "z", aspectratio = 1., size=(1280, 720))
# plot!.(camber_coords, color = :black, label = :none)
# [ mesh3d!(coord, colorscale = :viridis) for (coord, norm_Γ) in zip(horseshoe_coords, norm_Γs) ]
# plot!.(streams, color = :green, label = :none)

# gui();