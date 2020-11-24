### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ fb214f00-2e69-11eb-2128-ed0e241dc5b6
begin
	using Revise
	using StaticArrays
	using Rotations
	using BenchmarkTools
	using TimerOutputs
end

# ╔═╡ c37e7db0-2e6a-11eb-3f2b-116bc60bc1a0
begin
	include("D:/Academia/AeroMDAO/src/FoilParametrization.jl")
	using .FoilParametrization: naca4
	using AeroMDAO
end

# ╔═╡ 3a8ca800-2e75-11eb-3f19-13c43b57c698
using Plotly

# ╔═╡ fb446000-2e6c-11eb-1478-17896d7e8999
## Cuck

# ╔═╡ e4afd48e-2e72-11eb-2ca0-d799448b7e55
begin
	wing_foil = FoilParametrization.naca4((4,4,1,2))
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
	info(wing)
end

# ╔═╡ 73b21720-2e73-11eb-34e8-0309b1d0677a
begin
	htail_foil = FoilParametrization.naca4((0,0,1,2))
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
	info(htail)
end

# ╔═╡ e98c60e0-2e73-11eb-1180-7bdafb33f184
begin
	vtail_foil = FoilParametrization.naca4((0,0,0,9))
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
	info(vtail)
end

# ╔═╡ 2d7812e0-2e74-11eb-1a01-bbbcb4eeba21
wing_panels = paneller(wing, 10, 5);

# ╔═╡ 2f470680-2e74-11eb-2533-156fb8aebc2a
htail_panels = paneller(htail, 5, 5, translation = htail_location);

# ╔═╡ 2f47f0e2-2e74-11eb-2423-8573b83658ac
vtail1_panels = paneller(vtail, 5, 5, rotation = vtail1_angle, translation = vtail1_location);

# ╔═╡ 2f495070-2e74-11eb-3f6e-47dd24ce6362
vtail2_panels = paneller(vtail, 2, 2, rotation = vtail2_angle, translation = vtail2_location);

# ╔═╡ 382efe60-2e74-11eb-29f7-fb9dba6fdf9a
horseshoe_panels = [ wing_panels[1][:]; htail_panels[1][:]; vtail2_panels[1][:] ];

# ╔═╡ 40022ae0-2e74-11eb-1899-ef9068e78523
camber_panels = [ wing_panels[2][:]; htail_panels[2][:]; vtail2_panels[2][:] ];

# ╔═╡ 46ef3d1e-2e74-11eb-3860-6767f6de292f
ρ = 1.225

# ╔═╡ 492cd110-2e74-11eb-054d-61485a1035e8
ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.);

# ╔═╡ 5095ca12-2e74-11eb-079e-9370873531b0
Ω = SVector(0.0, 0.0, 0.0);

# ╔═╡ 55761c62-2e74-11eb-0347-ed05ff58c923
uniform = Freestream(10.0, 5.0, 0.0)

# ╔═╡ 5af96570-2e74-11eb-26f5-2b9f61ee57ae
@time force, drag, moment, horseshoes, Γs = solve_case(horseshoe_panels, camber_panels, uniform, Ω, ref, print = true) 

# ╔═╡ 174018a0-2e75-11eb-1904-8b3cb4197ccc
nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, Ω, uniform.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

# ╔═╡ 210418a0-2e75-11eb-16cd-f55588a20a5a
streams = plot_streamlines.(streamlines(uniform, Ω, horseshoe_panels, horseshoes, Γs, 5, 100));

# ╔═╡ 28a20450-2e75-11eb-0190-b9a22517983f
begin
	min_Γ, max_Γ = extrema(Γs)
	Γ_range = -map(-, min_Γ, max_Γ)
	norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ];
end

# ╔═╡ 310e3c80-2e75-11eb-3253-c16fa3c19458
begin
	aircraft_coords = plot_panels(horseshoe_panels)[:]
	camber_coords = plot_panels(camber_panels)[:]
	horseshoe_coords = plot_panels(horseshoe_panels);
	horseshoe_coords = plot_panels(vtail2_panels[1]);
end

# ╔═╡ 41469c50-2e75-11eb-2a31-1dd0d4048543
begin
	horse_xs = [ [ c[1] for c in panel ] for panel in horseshoe_coords ]
	horse_ys = [ [ c[2] for c in panel ] for panel in horseshoe_coords ]
	horse_zs = [ [ c[3] for c in panel ] for panel in horseshoe_coords ]

	camber_xs = [ [ c[1] for c in panel ] for panel in camber_coords ]
	camber_ys = [ [ c[2] for c in panel ] for panel in camber_coords ]
	camber_zs = [ [ c[3] for c in panel ] for panel in camber_coords ]

	streams_xs = [ [ c[1] for c in panel ] for panel in streams ]
	streams_ys = [ [ c[2] for c in panel ] for panel in streams ]
	streams_zs = [ [ c[3] for c in panel ] for panel in streams ];
end

# ╔═╡ 4ea7e752-2e75-11eb-13d0-356d2fabefbb
begin
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
end

# ╔═╡ 5f98e8c0-2e75-11eb-06e3-bf1763d368a7
Plotly.plot([ 
        [ trace for trace in trace_horses ]...,
        # [ trace for trace in trace_horsies ]..., 
        # [ trace for trace in trace_cambers ]...,
        # [ trace for trace in trace_streams ]...,
        # [ trace for trace in trace_aircraft ]...,
     ], 
     layout)

# ╔═╡ Cell order:
# ╠═fb446000-2e6c-11eb-1478-17896d7e8999
# ╠═fb214f00-2e69-11eb-2128-ed0e241dc5b6
# ╠═c37e7db0-2e6a-11eb-3f2b-116bc60bc1a0
# ╠═e4afd48e-2e72-11eb-2ca0-d799448b7e55
# ╠═73b21720-2e73-11eb-34e8-0309b1d0677a
# ╠═e98c60e0-2e73-11eb-1180-7bdafb33f184
# ╠═2d7812e0-2e74-11eb-1a01-bbbcb4eeba21
# ╠═2f470680-2e74-11eb-2533-156fb8aebc2a
# ╠═2f47f0e2-2e74-11eb-2423-8573b83658ac
# ╠═2f495070-2e74-11eb-3f6e-47dd24ce6362
# ╠═382efe60-2e74-11eb-29f7-fb9dba6fdf9a
# ╠═40022ae0-2e74-11eb-1899-ef9068e78523
# ╠═46ef3d1e-2e74-11eb-3860-6767f6de292f
# ╠═492cd110-2e74-11eb-054d-61485a1035e8
# ╠═5095ca12-2e74-11eb-079e-9370873531b0
# ╠═55761c62-2e74-11eb-0347-ed05ff58c923
# ╠═5af96570-2e74-11eb-26f5-2b9f61ee57ae
# ╠═174018a0-2e75-11eb-1904-8b3cb4197ccc
# ╠═210418a0-2e75-11eb-16cd-f55588a20a5a
# ╠═28a20450-2e75-11eb-0190-b9a22517983f
# ╠═310e3c80-2e75-11eb-3253-c16fa3c19458
# ╠═3a8ca800-2e75-11eb-3f19-13c43b57c698
# ╠═41469c50-2e75-11eb-2a31-1dd0d4048543
# ╠═4ea7e752-2e75-11eb-13d0-356d2fabefbb
# ╠═5f98e8c0-2e75-11eb-06e3-bf1763d368a7
