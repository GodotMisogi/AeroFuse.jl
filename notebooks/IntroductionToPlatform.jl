### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 4d7366f0-5ef1-11eb-2aa0-f981a375bc4d
using AeroMDAO

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
md"""Run `using AeroMDAO` in your REPL or in a Pluto notebook."""

# ╔═╡ fc63a290-5f07-11eb-04b0-8365bbe8d0db
md"Your first NACA-4 airfoil"

# ╔═╡ ba514670-5ef1-11eb-0dfa-75ddf5f9d1c7
foil = naca4((4,4,1,2));

# ╔═╡ f47f5f0e-5f07-11eb-0799-4b4ee6409a8c
begin 
	using Plots
	plot(first.(foil), last.(foil), aspectratio = 1)
end

# ╔═╡ f5b7b800-5f07-11eb-05c3-7bb00ffa7817
md"Let's plot it!"

# ╔═╡ 1d2e4200-5f08-11eb-3d11-b9b80171b08a
md"Here's how we define a `HalfWing`:"

# ╔═╡ 447d96f0-5ef2-11eb-235c-0b9c7abae1bd
md"""
```julia
HalfWing(foils  	:: Vector{Foil}, 
		 chords 	:: Vector{Real}, 
		 twists 	:: Vector{Real}, 
		 spans  	:: Vector{Real}, 
		 dihedrals  :: Vector{Real}, 
		 sweeps 	:: Vector{Real})
```
"""

# ╔═╡ 0c7e4130-5f08-11eb-2545-15359cc79c20
md"Here's how we define a `Vector` of `Foil`s:"

# ╔═╡ 0a1186f0-5f08-11eb-1350-571b8046b1b3
foils = Foil.(foil for i ∈ 1:3);

# ╔═╡ 51b9ea3e-5ef1-11eb-38e7-2563d31f4477
wing_right = HalfWing(foils,
                      [0.2, 0.2, 0.1],
					  [0., 2., 5.],
					  [1.705 / 2, 0.1],
					  [0., 60.],
					  [0., 30.]);

# ╔═╡ c7522a80-5f08-11eb-3a90-4b90f529f680
html"""<style>
main {
    max-width: 60%;
}
"""

# ╔═╡ Cell order:
# ╟─85fb2bf0-5eee-11eb-1726-b5c2d8f1a798
# ╟─9f4070c0-5eee-11eb-3e65-73511c496c89
# ╟─a90274f0-5eee-11eb-272f-9fea78ad929d
# ╟─11f8f630-5ef1-11eb-09df-9177abb61f44
# ╟─38ddc780-5ef1-11eb-22c8-792dbbe43c39
# ╠═4d7366f0-5ef1-11eb-2aa0-f981a375bc4d
# ╟─fc63a290-5f07-11eb-04b0-8365bbe8d0db
# ╠═ba514670-5ef1-11eb-0dfa-75ddf5f9d1c7
# ╟─f5b7b800-5f07-11eb-05c3-7bb00ffa7817
# ╠═f47f5f0e-5f07-11eb-0799-4b4ee6409a8c
# ╟─1d2e4200-5f08-11eb-3d11-b9b80171b08a
# ╟─447d96f0-5ef2-11eb-235c-0b9c7abae1bd
# ╟─0c7e4130-5f08-11eb-2545-15359cc79c20
# ╠═0a1186f0-5f08-11eb-1350-571b8046b1b3
# ╠═51b9ea3e-5ef1-11eb-38e7-2563d31f4477
# ╟─c7522a80-5f08-11eb-3a90-4b90f529f680
