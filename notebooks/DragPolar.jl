### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ ae5d8190-2f43-11eb-2c37-1d96d812053a
using LinearAlgebra

# ╔═╡ c695e3e0-2f3b-11eb-027f-ed0248d96007
using Plots

# ╔═╡ 70eded00-2f42-11eb-033d-b5c11b4da524
md"""
## Drag Polar
"""

# ╔═╡ 73f38150-2f41-11eb-0e87-d3b12942a42c
md"""
So this is how text rendering works, sheesh.
"""

# ╔═╡ b979b830-2f3b-11eb-2c1f-95035ea2a15b
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2

# ╔═╡ 7e002da0-2f42-11eb-33b8-edab9d45e9ee
md"""
And here's how math rendering works.
```math
\begin{equation*}
i\hbar = \frac{\partial \psi}{\partial t} \hat H
\end{equation*}
```
"""

# ╔═╡ bc8830b2-2f3b-11eb-1582-95999db29675
stall_speed(wing_loading, CL_max, ρ = 1.225) = 2 * wing_loading/(ρ * CL_max)

# ╔═╡ cb4cb52e-2f3b-11eb-0121-c92bc4b2b950
span_efficiency_factor(e, AR) = 1 / (π * e * AR)

# ╔═╡ d701fde0-2f3b-11eb-345f-955af341b009
drag_polar(CD0, k, CL, CL0 = 0.) = CD0 + k * (CL - CL0)^2

# ╔═╡ dca970c0-2f3b-11eb-21c6-9b2ceb7b62ca
CD0 = 1.5e-4

# ╔═╡ e19a16c0-2f3b-11eb-204b-a1d79f7a31ad
k = span_efficiency_factor(0.89, 6)

# ╔═╡ f5d84760-2f3b-11eb-3729-97a99b00087f
dragger(CL, CL0) = drag_polar(CD0, k, CL, CL0)

# ╔═╡ 106ffcd0-2f3c-11eb-21b8-c966a2a67288
dragger(0.5, 0.2)

# ╔═╡ 46d54fa0-2f3c-11eb-2c81-c972cc10fc41
cls = -1.5:0.05:1.5

# ╔═╡ 17dc7840-2f3c-11eb-26d4-67e3a57e9730
plot(dragger.(cls, -0.2), cls, marker=:dot)

# ╔═╡ 1849edd0-2f41-11eb-2d76-3f472da8b174
begin
    import Pkg
    Pkg.add(url="https://github.com/Pocket-titan/DarkMode.git")
    import DarkMode
	DarkMode.enable(theme="material-darker")
end

# ╔═╡ Cell order:
# ╟─70eded00-2f42-11eb-033d-b5c11b4da524
# ╠═ae5d8190-2f43-11eb-2c37-1d96d812053a
# ╟─73f38150-2f41-11eb-0e87-d3b12942a42c
# ╠═b979b830-2f3b-11eb-2c1f-95035ea2a15b
# ╟─7e002da0-2f42-11eb-33b8-edab9d45e9ee
# ╠═bc8830b2-2f3b-11eb-1582-95999db29675
# ╠═cb4cb52e-2f3b-11eb-0121-c92bc4b2b950
# ╠═d701fde0-2f3b-11eb-345f-955af341b009
# ╠═c695e3e0-2f3b-11eb-027f-ed0248d96007
# ╠═dca970c0-2f3b-11eb-21c6-9b2ceb7b62ca
# ╠═e19a16c0-2f3b-11eb-204b-a1d79f7a31ad
# ╠═f5d84760-2f3b-11eb-3729-97a99b00087f
# ╠═106ffcd0-2f3c-11eb-21b8-c966a2a67288
# ╠═46d54fa0-2f3c-11eb-2c81-c972cc10fc41
# ╠═17dc7840-2f3c-11eb-26d4-67e3a57e9730
# ╠═1849edd0-2f41-11eb-2d76-3f472da8b174
