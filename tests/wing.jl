## 
using Revise
includet("../src/MathTools.jl")
includet("../src/FoilParametrization.jl")

##
using StaticArrays
using .FoilParametrization: read_foil, foil_camthick, camthick_foil, cosine_foil, kulfan_CST, naca4
using .MathTools: linspace, tuparray, tupvector
using AeroMDAO
using DelimitedFiles
using Rotations

## Wing section setup
# alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
# alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
# alphas = [alpha_u alpha_l]
# dzs = (1e-4, 1e-4)
# foil = kulfan_CST(alphas, dzs, 0.2)

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
print_info(wing_right)

## Assembly
ρ = 1.225
ref = (0.25 * mean_aerodynamic_chord(wing_right), 0, 0)
uniform = Uniform(10.0, 5.0, 0.0)
@time horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing_right, uniform, ref, span_num = 10, chord_num = 5)

## Panel method: TO DO

wing_panels = mesh_wing(wing, 5, 5);

##
wing_coords = plot_panels(wing_panels)[:]
camber_coords = plot_panels(camber_panels)[:]
horseshoe_coords = plot_panels(horseshoe_panels)[:]
streams = plot_streamlines.(streamlines(uniform, horseshoe_panels, horseshoes, Γs, 3, 50));

##
min_Γ, max_Γ = extrema(Γs)
color_range = -map(-, min_Γ, max_Γ)
norm_Γs = [ 2(Γ - min_Γ)/color_range - 1 for Γ ∈ Γs ]

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


##
layout = Layout(
                title = "Penguins",
                scene=attr(aspectmode="manual", aspectratio=attr(x=1,y=1,z=1))
                )

trace_horses = [ mesh3d(
                        x = x,
                        y = y,
                        z = z,
                        intensity = repeat([norm_Γ], length(x)),
                        text = norm_Γ,
                        showscale = false,
                        ) for (x, y, z, norm_Γ) in zip(horse_xs, horse_ys, horse_zs, norm_Γs) ]

trace_horsies = [ scatter3d(
                            x = x,
                            y = y,
                            z = z,
                            mode = :lines, 
                            line = attr(color =:black),
                            showlegend = false,
                            ) for (x, y, z) in zip(horse_xs, horse_ys, horse_zs) ]

trace_cambers = [ scatter3d(
                       x = x,
                       y = y,
                       z = z,
                       mode = :lines, 
                       line = attr(color =:black),
                       showlegend = false,
                       ) for (x, y, z) in zip(camber_xs, camber_ys, camber_zs) ]

trace_streams = [ scatter3d(
                            x = x, 
                            y = y, 
                            z = z, 
                            mode = :lines, 
                            line = attr(color =:lightblue),
                            showlegend = false,
                            ) for (x, y, z) in zip(streams_xs, streams_ys, streams_zs) ]

plot([ 
        [ trace for trace in trace_horses ]...,
        [ trace for trace in trace_horsies ]..., 
        [ trace for trace in trace_cambers ]...,
        [ trace for trace in trace_streams ]... 
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