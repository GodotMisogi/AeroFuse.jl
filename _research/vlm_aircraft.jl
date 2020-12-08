## 
using Revise
using StaticArrays
using BenchmarkTools
using TimerOutputs
using ProfileView
using Rotations
using AeroMDAO

## Wing
wing_foil = naca4((4,4,1,2))
wing_right = HalfWing(Foil.(wing_foil for i ∈ 1:3), 
                      [0.18, 0.16, 0.08],
                      [0., 0., 0.],
                      [0.5, 0.5],
                      [0., 11.3],
                      [1.14, 8.])

wing = Wing(wing_right, wing_right)
println("Wing —")
print_info(wing)

## Horizontal tail
htail_foil = naca4((0,0,1,2))
htail_right = HalfWing(Foil.(htail_foil for i ∈ 1:2), 
                       [0.16, 0.08],
                       [0.0, 0.0],
                       [0.2],
                       [0.],
                       [30.])

htail = Wing(htail_right, htail_right)
println("Horizontal Tail —")
print_info(htail)

## Vertical tail
vtail_foil = naca4((0,0,0,9))
vtail = HalfWing(Foil.(vtail_foil for i ∈ 1:2), 
                 [0.08, 0.02],
                 [0.0, 0.0],
                 [0.04],
                 [0.],
                 [60.])
println("Vertical Tail —")
print_info(vtail)

## Panelling
htail_location = [1., 0., 0.]

vtail1_angle = AngleAxis{Float64}(π/4, 1, 0, 0)
vtail1_location = [1 + 0.2 * tan(π/6), 0.2, 0]

vtail2_angle = AngleAxis{Float64}(3π/4, 1, 0, 0)
vtail2_location = [1 + 0.2 * tan(π/6), -0.2, 0]

wing_panels = paneller(wing, 10, 5)
htail_panels = paneller(htail, 5, 5, translation = htail_location)
vtail1_panels = paneller(vtail, 2, 2, rotation = vtail1_angle, translation = vtail1_location)
vtail2_panels = paneller(vtail, 2, 2, rotation = vtail2_angle, translation = vtail2_location);

##
horseshoe_panels = [ wing_panels[1][:]; htail_panels[1][:]; vtail1_panels[1][:]; vtail2_panels[1][:] ]
camber_panels = [ wing_panels[2][:]; htail_panels[2][:]; vtail1_panels[2][:]; vtail2_panels[2][:] ];

## Assembly
reset_timer!()

ρ = 1.225
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
Ω = SVector(0.0, 0.0, 0.0)
freestream = Freestream(10.0, 0.0, 0.0, Ω)
@time force, drag, moment, horseshoes, Γs = solve_case(horseshoe_panels, camber_panels, freestream, ref, print = true) 

@timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, freestream.Ω, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

println("Nearfield:") 
print_dynamics(nearfield_coeffs...)

print_timer();

## Streamlines
reset_timer!()

@timeit "Computing Streamlines" streams = plot_streamlines.(streamlines(freestream, Ω, horseshoe_panels, horseshoes, Γs, 5, 100));

print_timer()

## Normalisation
min_Γ, max_Γ = extrema(Γs)
Γ_range = -map(-, min_Γ, max_Γ)
norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ]

##
# aircraft_coords = plot_panels(horseshoe_panels)[:]
camber_coords = plot_panels(camber_panels)[:]
horseshoe_coords = plot_panels(horseshoe_panels)[:];
# horseshoe_coords = plot_panels(vtail2_panels[1]);

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
        # [ trace for trace in trace_horsies ]..., 
        # [ trace for trace in trace_cambers ]...,
        [ trace for trace in trace_streams ]...,
        # [ trace for trace in trace_aircraft ]...,
        [ trace for trace in trace_horses ]...,
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