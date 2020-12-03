## 
using Revise
using StaticArrays
using BenchmarkTools
using TimerOutputs
using ProfileView

## Custom packages
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")
using .FoilParametrization: read_foil, kulfan_CST, naca4
using AeroMDAO
using Rotations

## Wing section setup
# alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
# alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
# alphas = [alpha_u alpha_l]
# dzs = (1e-4, 1e-4)
# foil = kulfan_CST(alphas, dzs, 0.2)

# Wing setup
wing_foil = naca4((4,4,1,2))
wing_num_secs = 3
wing_foils = [ wing_foil for i ∈ 1:wing_num_secs ]

airfoils = Foil.(wing_foils)
wing_chords = [0.18, 0.16, 0.08]
wing_twists = [0., 0., 0.]
wing_spans = [0.5, 0.5]
wing_dihedrals = [0., 11.3]
wing_sweeps = [1.14, 8.]

wing_right = HalfWing(airfoils, wing_chords, wing_spans, wing_dihedrals, wing_sweeps, wing_twists)
wing = Wing(wing_right, wing_right)
println("Wing —")
print_info(wing)

## Horizontal tail
htail_foil = naca4((0,0,1,2))
htail_num_secs = 2
htail_foils = [ htail_foil for i ∈ 1:htail_num_secs ]

htail_airfoils = Foil.(htail_foils)
htail_chords = [0.16, 0.08]
htail_twists = [0.0, 0.0]
htail_spans = [0.2]
htail_dihedrals = [0]
htail_sweeps = [30]
htail_location = [1,0,0]

htail_right = HalfWing(htail_airfoils, htail_chords, htail_spans, htail_dihedrals, htail_sweeps, htail_twists)
htail = Wing(htail_right, htail_right)
println("Horizontal Tail —")
print_info(htail)

## Vertical tail
vtail_foil = naca4((0,0,0,9))
vtail_num_secs = 2
vtail_foils = [ vtail_foil for i ∈ 1:vtail_num_secs ]

vtail_airfoils = Foil.(vtail_foils)
vtail_chords = [0.08, 0.02]
vtail_twists = [0.0, 0.0]
vtail_spans = [0.04]
vtail_dihedrals = [0.0]
vtail_sweeps = [60]

vtail1_angle = AngleAxis{Float64}(π/2, 1, 0, 0)
vtail1_location = [1 + 0.2 * tan(π/6), 0.2, 0]

vtail2_angle = AngleAxis{Float64}(π/2, 1, 0, 0)
vtail2_location = [1 + 0.2 * tan(π/6), -0.2, 0]

vtail = HalfWing(vtail_airfoils, vtail_chords, vtail_spans, vtail_dihedrals, vtail_sweeps, vtail_twists)
println("Vertical Tail —")
print_info(vtail)

## Panelling
wing_panels = paneller(wing, 10, 5)
htail_panels = paneller(htail, 5, 5, translation = htail_location)
vtail1_panels = paneller(vtail, 5, 5, rotation = vtail1_angle, translation = vtail1_location)
vtail2_panels = paneller(vtail, 2, 2, rotation = vtail2_angle, translation = vtail2_location);

##
horseshoe_panels = [ wing_panels[1][:]; htail_panels[1][:]; vtail2_panels[1][:] ]
camber_panels = [ wing_panels[2][:]; htail_panels[2][:]; vtail2_panels[2][:] ];

## Assembly
reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
uniform = Freestream(10.0, 5.0, 0.0)
@time force, drag, moment, horseshoes, Γs = solve_case(horseshoe_panels, camber_panels, uniform, Ω, ref, print = true) 

@timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, Ω, uniform.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

println("Nearfield:") 
print_dynamics(nearfield_coeffs...)

print_timer();
## Streamlines
reset_timer!()

@timeit "Computing Streamlines" streams = plot_streamlines.(streamlines(uniform, Ω, horseshoe_panels, horseshoes, Γs, 5, 100));

print_timer()

##
min_Γ, max_Γ = extrema(Γs)
Γ_range = -map(-, min_Γ, max_Γ)
norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ];

aircraft_coords = plot_panels(horseshoe_panels)[:]
camber_coords = plot_panels(camber_panels)[:]
horseshoe_coords = plot_panels(horseshoe_panels);
horseshoe_coords = plot_panels(vtail2_panels[1]);

##
using PlotlyJS

##
horse_xs = [ [ c[1] for c in panel ] for panel in horseshoe_coords ]
horse_ys = [ [ c[2] for c in panel ] for panel in horseshoe_coords ]
horse_zs = [ [ c[3] for c in panel ] for panel in horseshoe_coords ]

camber_xs = [ [ c[1] for c in panel ] for panel in camber_coords ]
camber_ys = [ [ c[2] for c in panel ] for panel in camber_coords ]
camber_zs = [ [ c[3] for c in panel ] for panel in camber_coords ]

streams_xs = [ [ c[1] for c in panel ] for panel in streams ]
streams_ys = [ [ c[2] for c in panel ] for panel in streams ]
streams_zs = [ [ c[3] for c in panel ] for panel in streams ];

# aircraft_xs = [ [ c[1] for c in panel ] for panel in aircraft_coords ]
# aircraft_ys = [ [ c[2] for c in panel ] for panel in aircraft_coords ]
# aircraft_zs = [ [ c[3] for c in panel ] for panel in aircraft_coords ];

##
layout = Layout(
                title = "Penguins",
                scene=attr(aspectmode="manual", aspectratio=attr(x=1,y=1,z=1)),
                zlim=(-0.1, 5.0)
                )

trace_horses = [ PlotlyJS.mesh3d(
                        x = x,
                        y = y,
                        z = z,
                        intensity = repeat([norm_Γ], length(x)),
                        text = norm_Γ,
                        showscale = false,
                        ) for (x, y, z, norm_Γ) in zip(horse_xs, horse_ys, horse_zs, norm_Γs) ]

trace_horsies = [ PlotlyJS.scatter3d(
                            x = x,
                            y = y,
                            z = z,
                            mode = :lines, 
                            line = attr(color = :black),
                            showlegend = false,
                            ) for (x, y, z) in zip(horse_xs, horse_ys, horse_zs) ]

trace_cambers = [ PlotlyJS.scatter3d(
                       x = x,
                       y = y,
                       z = z,
                       mode = :lines, 
                       line = attr(color = :black),
                       showlegend = false,
                       ) for (x, y, z) in zip(camber_xs, camber_ys, camber_zs) ]

trace_streams = [ PlotlyJS.scatter3d(
                            x = x, 
                            y = y, 
                            z = z, 
                            mode = :lines, 
                            line = attr(color = :lightblue),
                            showlegend = false,
                            ) for (x, y, z) in zip(streams_xs, streams_ys, streams_zs) ];

# trace_aircraft =    [ PlotlyJS.scatter3d(
#                             x = x,
#                             y = y,
#                             z = z,
#                             mode = :lines, 
#                             line = attr(color = :black),
#                             showlegend = false,
#                             ) for (x, y, z) in zip(aircraft_xs, aircraft_ys, aircraft_zs) ];

PlotlyJS.plot([ 
        [ trace for trace in trace_horses ]...,
        # [ trace for trace in trace_horsies ]..., 
        # [ trace for trace in trace_cambers ]...,
        # [ trace for trace in trace_streams ]...,
        # [ trace for trace in trace_aircraft ]...,
     ], 
     layout)


## Plotting
# using Plots
# plotlyjs()

# ##

# plot(xaxis = "x", yaxis = "y", zaxis = "z", aspectratio = 1., zlim=(-0.1, 3.0), size=(1280, 720))
# plot!.(camber_coords, color = :black, label = :none)
# plot!.(horseshoe_coords, color = :blue, label = :none)
# plot!.(streams, color = :green, label = :none)

# gui();