### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f7132010-404c-11eb-157e-0dfc3c36ead3
begin
	using PlutoUI
	using Revise
	using StaticArrays
	using AeroMDAO
	using Plots
	gr()
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

# ╔═╡ 0ac53d50-404d-11eb-0485-637991172933
begin
	ρ = 1.225
	ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
	Ω = SVector(0.0, 0.0, 0.0)
	uniform = Freestream(10.0, 5.0, 5.0, Ω)
end;

# ╔═╡ 0add0b10-404d-11eb-0b56-91f9e3431c70
@time coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 10, print = true);

# ╔═╡ 81833050-404d-11eb-0dbd-7fc5e4d69d77
(collect ∘ zip)(["CL", "CDi", "CY", "Cl", "Cm", "Cn", "p", "q", "r"], coeffs)

# ╔═╡ 0b11fdc2-404d-11eb-269c-6904c57e599d
begin
	horseshoe_coords = plot_panels(horseshoe_panels[:])
	# horses_coords = plot_panels(horseshoe_panels[:], Γs[:])
	camber_coords = plot_panels(camber_panels[:])
	wing_coords = plot_surface(wing);
end;

# ╔═╡ 158416d0-40f7-11eb-34a9-c193e6456ff3
begin
	φ_s = @bind φ Slider(0:1e-2:90, default = 15)
	ψ_s = @bind ψ Slider(0:1e-2:90, default = 30)
	z_s = @bind z_limit Slider(0:1e-2:3, default = 1.)
	md"""
	Horizontal: $(φ_s)
	Vertical: $(ψ_s)
	
	*z*-scale: $(z_s)
	"""
end

# ╔═╡ 0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
begin
	plot(xaxis = "x", yaxis = "y", zaxis = "z",
		 aspect_ratio = 1, 
		 camera = (φ, ψ),
		 zlim = (-0.1, z_limit))
	plot!.(horseshoe_coords, color = :black, label = :none)
	# plot!.(streams, color = :green, label = :none)
	plot!()
end

# ╔═╡ Cell order:
# ╠═f7132010-404c-11eb-157e-0dfc3c36ead3
# ╠═0a99e790-404d-11eb-2b93-f3d0acd94967
# ╠═0aa950e0-404d-11eb-3a8a-bfe40df19222
# ╠═0aae0bd0-404d-11eb-1245-9186c5b9c21f
# ╠═0aaf4450-404d-11eb-31ba-757f8422c1ab
# ╠═0ac53d50-404d-11eb-0485-637991172933
# ╠═0add0b10-404d-11eb-0b56-91f9e3431c70
# ╠═81833050-404d-11eb-0dbd-7fc5e4d69d77
# ╠═0b11fdc2-404d-11eb-269c-6904c57e599d
# ╟─158416d0-40f7-11eb-34a9-c193e6456ff3
# ╟─0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
