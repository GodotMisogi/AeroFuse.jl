### A Pluto.jl notebook ###
# v0.12.18

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

# ╔═╡ 4d7366f0-5ef1-11eb-2aa0-f981a375bc4d
using AeroMDAO

# ╔═╡ 92bbf44e-5fb6-11eb-1c27-6ba435bca7ba
using PlutoUI

# ╔═╡ 85fb2bf0-5eee-11eb-1726-b5c2d8f1a798
md"# MECH 3620 -- Introduction to Aircraft Design Platform"

# ╔═╡ 9f4070c0-5eee-11eb-3e65-73511c496c89
md"""## Fixed Wing-VTOL Hybrid Case

![](https://godot-bloggy.xyz/DesignFramework.svg)

"""

# ╔═╡ a90274f0-5eee-11eb-272f-9fea78ad929d
md"""
**References**:

1. Tyan, Maxim et al. _Comprehensive preliminary sizing/resizing method for a fixed wing – VTOL electric UAV_. Aerospace Science and Technology, 2017.
2. Drela, Mark. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
"""

# ╔═╡ 11f8f630-5ef1-11eb-09df-9177abb61f44
md"""
## Your First Wing
"""

# ╔═╡ 38ddc780-5ef1-11eb-22c8-792dbbe43c39
md"""Run the following line in your REPL or in a Pluto notebook."""

# ╔═╡ fc63a290-5f07-11eb-04b0-8365bbe8d0db
md"Your first NACA-4 airfoil"

# ╔═╡ 3b18703e-5f28-11eb-2cf2-21542c05c32a
digits = (2,4,1,2);

# ╔═╡ ba514670-5ef1-11eb-0dfa-75ddf5f9d1c7
airfoil = naca4(digits);

# ╔═╡ f47f5f0e-5f07-11eb-0799-4b4ee6409a8c
begin 
	using Plots
	plot(getindex.(airfoil, 1), getindex.(airfoil, 2), marker=2, aspectratio = 1, label = "NACA $(digits...)")
end

# ╔═╡ f5b7b800-5f07-11eb-05c3-7bb00ffa7817
md"Let's plot it!"

# ╔═╡ 1d2e4200-5f08-11eb-3d11-b9b80171b08a
md"Here's how we define a `HalfWing` in `AeroMDAO`:"

# ╔═╡ 447d96f0-5ef2-11eb-235c-0b9c7abae1bd
md"""
```julia
HalfWing(foils  	:: Vector{Foil}, 	# Foil profiles
		 chords 	:: Vector{Real}, 	# Chord lengths
		 twists 	:: Vector{Real}, 	# Twist angles
		 spans  	:: Vector{Real}, 	# Section span lengths
		 dihedrals  :: Vector{Real}, 	# Dihedral angles
		 sweeps 	:: Vector{Real})	# Leading-edge sweep angles
```
"""

# ╔═╡ 0c7e4130-5f08-11eb-2545-15359cc79c20
md"Here's how we define a `Vector` of `Foil`s. Note the use of broadcasting and comprehensions."

# ╔═╡ 0a1186f0-5f08-11eb-1350-571b8046b1b3
foils = Foil.(airfoil for i ∈ 1:3)

# ╔═╡ 45c38130-5fb2-11eb-1c9c-4367ef32d70c
md"Now that we have our foil profiles in the appropriate type, we can define our wing:"

# ╔═╡ 51b9ea3e-5ef1-11eb-38e7-2563d31f4477
wing_right = HalfWing(foils,
                      [0.4, 0.2, 0.1],
					  [0., 2., 5.],
					  [1.0, 0.1],
					  [0., 60.],
					  [0., 30.]);

# ╔═╡ ce8598f0-5f26-11eb-12b4-95bcf7d3e67e
md"We can create a `Wing` by feeding two `HalfWing`s to it:
```julia
Wing(left :: HalfWing, right :: HalfWing)
```
"

# ╔═╡ fe187750-5f2f-11eb-023a-f9d013d930b7
md"In this case, we'd like a symmetric wing, so we just feed `wing_right` as both arguments."

# ╔═╡ f470a0f0-5f26-11eb-32ed-97bff0246bfd
wing = Wing(wing_right, wing_right);

# ╔═╡ effac010-5fb1-11eb-1262-0b623dac72f0
md"Here's what the outline of our wing looks like:"

# ╔═╡ 809495c2-5fb6-11eb-2b5c-05f5c784c45f
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

# ╔═╡ 984eebb0-5f30-11eb-17e6-21280828a17a
begin
	cuck = plot(xaxis = "x", yaxis = "y", zaxis = "z",
		 aspect_ratio = 1, 
		 camera = (φ, ψ),
		 zlim = (-0.1, z_limit)
		)
	wing_coords = plot_wing(wing)
	plot!(cuck, wing_coords, label = "Wing")
end

# ╔═╡ c7522a80-5f08-11eb-3a90-4b90f529f680
html"""<style>
main {
    max-width: 60%;
}
"""

# ╔═╡ 04e44e50-5f22-11eb-1495-316729d2795a
begin
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));
	definition(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("warning", "Definition", [text]));
	exercise(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("correct", "Exercise", [text]));
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
end;

# ╔═╡ 2c9a11f0-5f27-11eb-2932-259cf298a539
definition(md"A _wing section_ consists of two foil profiles and their chord lengths and twist angles. Between them is their span length with associated dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for `spans`, `dihedrals`, `sweeps`, and ``n+1`` entries for `foils`, `chords`, and `twists` for some ``n \in \mathbb Z``.")

# ╔═╡ ae5493b0-5f21-11eb-133b-7de74aefd6fa
exercise(md"Create a `Wing` consisting of 4 spanwise sections with different NACA 4-digit foil profiles using whatever dimensions you like, as long as they're physically reasonable.")

# ╔═╡ Cell order:
# ╟─85fb2bf0-5eee-11eb-1726-b5c2d8f1a798
# ╟─9f4070c0-5eee-11eb-3e65-73511c496c89
# ╟─a90274f0-5eee-11eb-272f-9fea78ad929d
# ╟─11f8f630-5ef1-11eb-09df-9177abb61f44
# ╟─38ddc780-5ef1-11eb-22c8-792dbbe43c39
# ╠═4d7366f0-5ef1-11eb-2aa0-f981a375bc4d
# ╟─fc63a290-5f07-11eb-04b0-8365bbe8d0db
# ╠═3b18703e-5f28-11eb-2cf2-21542c05c32a
# ╠═ba514670-5ef1-11eb-0dfa-75ddf5f9d1c7
# ╟─f5b7b800-5f07-11eb-05c3-7bb00ffa7817
# ╟─f47f5f0e-5f07-11eb-0799-4b4ee6409a8c
# ╟─1d2e4200-5f08-11eb-3d11-b9b80171b08a
# ╟─2c9a11f0-5f27-11eb-2932-259cf298a539
# ╟─447d96f0-5ef2-11eb-235c-0b9c7abae1bd
# ╟─0c7e4130-5f08-11eb-2545-15359cc79c20
# ╠═0a1186f0-5f08-11eb-1350-571b8046b1b3
# ╟─45c38130-5fb2-11eb-1c9c-4367ef32d70c
# ╠═51b9ea3e-5ef1-11eb-38e7-2563d31f4477
# ╟─ce8598f0-5f26-11eb-12b4-95bcf7d3e67e
# ╟─fe187750-5f2f-11eb-023a-f9d013d930b7
# ╠═f470a0f0-5f26-11eb-32ed-97bff0246bfd
# ╟─effac010-5fb1-11eb-1262-0b623dac72f0
# ╟─809495c2-5fb6-11eb-2b5c-05f5c784c45f
# ╟─984eebb0-5f30-11eb-17e6-21280828a17a
# ╟─ae5493b0-5f21-11eb-133b-7de74aefd6fa
# ╟─c7522a80-5f08-11eb-3a90-4b90f529f680
# ╟─04e44e50-5f22-11eb-1495-316729d2795a
# ╟─92bbf44e-5fb6-11eb-1c27-6ba435bca7ba
