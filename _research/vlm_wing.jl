## 
using BenchmarkTools
# using ProfileView
using Revise
using StaticArrays
using TimerOutputs
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
reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 5.0, 0.0, Ω)
@time nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10) 

print_timer();

begin
    println("\nNearfield:")
    print_dynamics(nf_coeffs...)
    println("\nFarfield:")
    print_dynamics(ff_coeffs...)
end

##
using PlotlyJS

##
layout = Layout(title = "Vortex Lattice",
                scene = attr(aspectratio=attr(x=1,y=1,z=1)));

## Streamlines
reset_timer!()

num_points = 30
max_z = 0.02
y = span(wing) / 2 - 0.5
# seed = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))

span_points = 10
init = trailing_chopper(wing.right, span_points) 
dx, dy, dz = 0, 0, 1e-3
seed = [ init .+ Ref([dx, dy, dz]); 
         init .+ Ref([dx, dy, -dz]) ]

@timeit "Computing Streamlines" trace_streams = trace_streamlines(uniform, seed, horseshoes[:], Γs[:], 2, 100);

print_timer()

##
trace_horsies = trace_panels(horseshoe_panels[:])
trace_horses = trace_panels(horseshoe_panels[:], Γs[:])
trace_cambers = trace_panels(camber_panels[:])
trace_wing = trace_surface(wing)

PlotlyJS.plot(
            [
                [ trace for trace in trace_horsies ]...,
                [ trace for trace in trace_horses ]...,
                [ trace for trace in trace_streams ]...,
                [ trace for trace in trace_cambers ]...,
                # [ trace for trace in trace_wing ]...,
            ],
            layout)