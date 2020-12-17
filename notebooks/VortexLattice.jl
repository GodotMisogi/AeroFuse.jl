### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ f7132010-404c-11eb-157e-0dfc3c36ead3
begin
	using Revise
	using StaticArrays
	using AeroMDAO
	using Plots
	plotlyjs()
end;

# ╔═╡ 0a99e790-404d-11eb-2b93-f3d0acd94967
foil = naca4((4,4,1,2));

# ╔═╡ 0aa950e0-404d-11eb-3a8a-bfe40df19222
wing_right = HalfWing(Foil.(foil for i ∈ 1:3),    # Foils
                      [0.18, 0.16, 0.08],         # Chords
                      [2., 0., -2.],              # Twists
                      [0.5, 0.2],                 # Spans
                      [0., 11.3],                 # Dihedrals
                      [1.14, 8.]);                # Sweeps

# ╔═╡ 0aae0bd0-404d-11eb-1245-9186c5b9c21f
wing = Wing(wing_right, wing_right);

# ╔═╡ 0aaf4450-404d-11eb-31ba-757f8422c1ab
info(wing)

# ╔═╡ c8d66030-4061-11eb-1ce5-a3fef1793d0d
wing_coords = plot_surface(wing);

# ╔═╡ 0ac53d50-404d-11eb-0485-637991172933
begin
	ρ = 1.225
	ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
	Ω = SVector(0.0, 0.0, 0.0)
	uniform = Freestream(10.0, 5.0, 0.0, Ω)
end;

# ╔═╡ 0add0b10-404d-11eb-0b56-91f9e3431c70
@time coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10, print = true);

# ╔═╡ 81833050-404d-11eb-0dbd-7fc5e4d69d77
coeffs

# ╔═╡ 0b11fdc2-404d-11eb-269c-6904c57e599d
begin
	horseshoe_coords = plot_panels(horseshoe_panels[:])
	# horses_coords = plot_panels(horseshoe_panels[:], Γs[:])
	camber_coords = plot_panels(camber_panels[:])
end;

# ╔═╡ 0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
begin
	plot(xaxis = "x", yaxis = "y", zaxis = "z",
		 aspect_ratio = 1, 
		 camera = (30, 30))
	plot!.(horseshoe_coords, color = :black, label = :none)
	# plot!.(streams, color = :green, label = :none)
	plot!()
end

# ╔═╡ 532dc2e0-4063-11eb-1405-afd6744a5c61
streams = streamlines(uniform, horseshoe_panels[:], horseshoes[:], Γs[:], 2, 100);

# ╔═╡ Cell order:
# ╠═f7132010-404c-11eb-157e-0dfc3c36ead3
# ╠═0a99e790-404d-11eb-2b93-f3d0acd94967
# ╠═0aa950e0-404d-11eb-3a8a-bfe40df19222
# ╠═0aae0bd0-404d-11eb-1245-9186c5b9c21f
# ╠═0aaf4450-404d-11eb-31ba-757f8422c1ab
# ╠═c8d66030-4061-11eb-1ce5-a3fef1793d0d
# ╠═0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
# ╠═0ac53d50-404d-11eb-0485-637991172933
# ╠═0add0b10-404d-11eb-0b56-91f9e3431c70
# ╠═81833050-404d-11eb-0dbd-7fc5e4d69d77
# ╠═0b11fdc2-404d-11eb-269c-6904c57e599d
# ╠═532dc2e0-4063-11eb-1405-afd6744a5c61
