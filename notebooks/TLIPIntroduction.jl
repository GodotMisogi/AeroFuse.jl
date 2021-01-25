### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 25d361d0-407a-11eb-2e31-e393bd8aa7e1
md"""
# TLIP - Aircraft Design
"""

# ╔═╡ 670c4350-408f-11eb-0ded-f7d2115bd86c
md"""### Fixed Wing-VTOL Hybrid Case

Example Case:
![](https://godot-bloggy.xyz/DesignFramework.svg)

"""

# ╔═╡ 0700cbf0-40f5-11eb-1243-47e213dd19fb
md"""
**References**:

1. Tyan, Maxim et al. _Comprehensive preliminary sizing/resizing method for a fixed wing – VTOL electric UAV_. Aerospace Science and Technology, 2017.
2. Drela, Mark. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
"""

# ╔═╡ 2a2ab690-4091-11eb-0c34-89a03ac496c5
md"""## Progress

* Development of notebooks for experiments in MECH 3620 -- Aircraft Design. ✓
* Development of initial weight estimation code. ✓
* Development of constraint analysis code. ✓
* Development of tail sizing code. ✓
* Development of doublet-source panel method code for 2D aerodynamics analyses. ✓
* Development of vortex-lattice method code for 3D aerodynamics analyses. ✓
"""

# ╔═╡ 067fe520-4092-11eb-2cb9-c70f29ceb3da
md"""## To Do

1. Formulate a 'simple' optimization problem for MECH 3620.
2. Evaluate derivatives for dynamic stability using vortex-lattice code for MECH 3670.
3. Develop code for performance evaluations of aircraft for MECH 3670. 
4. Develop a module for structural analyses of aircraft. 
"""

# ╔═╡ 00a45890-4103-11eb-2fcf-cdb152320810


# ╔═╡ fbcffce2-40f1-11eb-2eee-87b8f493931c
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ edc6d920-40f1-11eb-04af-cf216294003a
hint(md"""
	It can also use HTML like this!
	""")

# ╔═╡ Cell order:
# ╠═25d361d0-407a-11eb-2e31-e393bd8aa7e1
# ╟─edc6d920-40f1-11eb-04af-cf216294003a
# ╟─670c4350-408f-11eb-0ded-f7d2115bd86c
# ╟─0700cbf0-40f5-11eb-1243-47e213dd19fb
# ╟─2a2ab690-4091-11eb-0c34-89a03ac496c5
# ╟─067fe520-4092-11eb-2cb9-c70f29ceb3da
# ╟─00a45890-4103-11eb-2fcf-cdb152320810
# ╟─fbcffce2-40f1-11eb-2eee-87b8f493931c
