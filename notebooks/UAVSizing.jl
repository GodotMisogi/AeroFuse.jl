### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# ╔═╡ 23137140-3091-11eb-3f8a-2b8bd4b426b5
md"""
## UAV Sizing
"""

# ╔═╡ 5fb5d550-3093-11eb-0691-09050e6939eb
md"""
```math
\begin{equation}
S_{HT} = \frac{C_{HT} \bar{c}_w S_w}{l_{HT}}
\end{equation}
```
"""

# ╔═╡ 304ebdb0-3091-11eb-0603-dd2131ef1226
horizontal_tail_area(V_H, wing_mac, wing_area, l_h) = V_H * wing_mac * wing_area / l_h

# ╔═╡ f5ae15b2-3091-11eb-3c08-53653717475a
rect_chord(area, span) = area / span

# ╔═╡ 2eacd810-3092-11eb-2cbb-d78bea2a7245
vertical_tail_area(V_V, wing_span, wing_area, l_v) = V_V * wing_span * wing_area / 2l_v

# ╔═╡ 40f52c20-3092-11eb-19c9-a3c63c15173b
vertical_tail_span(vtail_area, htail_chord, η_VT) = 2 * vtail_area / htail_chord / (1 + 1 / η_VT)

# ╔═╡ 95653b10-3092-11eb-2544-f1ec0d7148f3
htail_arm(x_VT, vtail_span, vtail_LE_sweep, htail_chord, x_CG) = x_VT + vtail_span * vtail_LE_sweep + 0.25 * htail_chord - x_CG

# ╔═╡ c1d62f10-3092-11eb-3421-355fe114bda5
vtail_arm(x_MAC_VT, vtail_mac) = x_MAC_VT+ 0.25 * vtail_mac

# ╔═╡ d5b5ec00-3092-11eb-2d10-87fcea8784a9
x_mac(x_CG, vtail_area, vtail_root, vtail_tip) = x_CG + vtail_area * (vtail_root + 2 * vtail_tip) / 3(vtail_root + vtail_tip)

# ╔═╡ 03ad8000-3093-11eb-1068-19d569a57d29
begin
	x_CG = 1.
	x_VT = 3.
	l_VT = x_VT - x_CG
end

# ╔═╡ 5c3009f0-3093-11eb-0ce2-9de2e3be4541


# ╔═╡ 18c1fe50-3091-11eb-060c-c944f6e2e2c3
begin
    import Pkg
    Pkg.add(url="https://github.com/Pocket-titan/DarkMode.git")
    import DarkMode
	DarkMode.enable(theme="material-darker")
end

# ╔═╡ Cell order:
# ╟─23137140-3091-11eb-3f8a-2b8bd4b426b5
# ╟─5fb5d550-3093-11eb-0691-09050e6939eb
# ╠═304ebdb0-3091-11eb-0603-dd2131ef1226
# ╠═f5ae15b2-3091-11eb-3c08-53653717475a
# ╠═2eacd810-3092-11eb-2cbb-d78bea2a7245
# ╠═40f52c20-3092-11eb-19c9-a3c63c15173b
# ╠═95653b10-3092-11eb-2544-f1ec0d7148f3
# ╠═c1d62f10-3092-11eb-3421-355fe114bda5
# ╠═d5b5ec00-3092-11eb-2d10-87fcea8784a9
# ╠═03ad8000-3093-11eb-1068-19d569a57d29
# ╠═5c3009f0-3093-11eb-0ce2-9de2e3be4541
# ╟─18c1fe50-3091-11eb-060c-c944f6e2e2c3
