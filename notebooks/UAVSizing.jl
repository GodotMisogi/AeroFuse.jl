### A Pluto.jl notebook ###
# v0.12.14

using Markdown
using InteractiveUtils

# ╔═╡ 23137140-3091-11eb-3f8a-2b8bd4b426b5
md"""
## UAV Sizing
"""

# ╔═╡ cbac90c0-3262-11eb-1908-87d130f6d01c
md"""
Tail Sizing
"""

# ╔═╡ 5fb5d550-3093-11eb-0691-09050e6939eb
md"""
```math
\begin{equation}
S_{HT} = \frac{C_{HT} \bar{c}_w S_w}{l_{HT}}
\end{equation}
```
"""

# ╔═╡ df735200-3263-11eb-3c15-c9f42ce9fbef
md"""
Tail sizing functions:
"""

# ╔═╡ 304ebdb0-3091-11eb-0603-dd2131ef1226
begin
	horizontal_tail_area(V_H, wing_mac, wing_area, l_h) = V_H * wing_mac * wing_area / l_h
	vertical_tail_area(V_V, wing_span, wing_area, l_v) = V_V * wing_span * wing_area / 2l_v
	vertical_tail_span(vtail_area, htail_chord, η_VT) = 2 * vtail_area / htail_chord / (1 + 1 / η_VT)
	htail_arm(x_VT, vtail_span, vtail_LE_sweep, htail_chord, x_CG) = x_VT + vtail_span * vtail_LE_sweep + 0.25 * htail_chord - x_CG
	vtail_arm(x_MAC_VT, vtail_mac) = x_MAC_VT+ 0.25 * vtail_mac
end

# ╔═╡ f5ae15b2-3091-11eb-3c08-53653717475a
begin
	rect_chord(area, span) = area / span
	mean_aerodynamic_chord(c_r, λ) = 2/3 * c_r * (1 + λ + λ^2) / (1 + λ)
	aspect_ratio(span, area) = span^2 / area
end

# ╔═╡ d5b5ec00-3092-11eb-2d10-87fcea8784a9
x_mac(x_CG, vtail_area, vtail_root, vtail_tip) = x_CG + vtail_area * (vtail_root + 2 * vtail_tip) / 3(vtail_root + vtail_tip)

# ╔═╡ 8cfe52f0-3267-11eb-2b11-23458f7f2966
md"""
Wing Parameters
"""

# ╔═╡ 7f48f152-3263-11eb-36d1-cfd2a24b2cd9
begin
	x_w = 0.33
	Λ_LE_w = deg2rad(0)
	Λ_TE_w = deg2rad(0)
	c_r_w = 0.20
	λ_w = 1.0
	c_t_w = λ_w * c_r_w
	mac_w = mean_aerodynamic_chord(c_r_w, λ_w)
	b_w = 1.0
	S_w = (1 + λ_w) / 2 * c_r_w
end

# ╔═╡ 98a56ace-3267-11eb-3627-ed82836734c2
md"""
Propeller locations
"""

# ╔═╡ 3d4d2870-3263-11eb-3e80-1f3f5db82243
begin
	D_VTOL_prop = 0.1
	D_FW_prop = 0.1
	c_prop = 0.01
end

# ╔═╡ 39d637e0-3268-11eb-2536-033d1b0d580e
md"""
Horizontal tail.
"""

# ╔═╡ f11e9b60-3266-11eb-0c1f-fd0eaea6cc5d
V_H = 0.5

# ╔═╡ 5c3009f0-3093-11eb-0ce2-9de2e3be4541
begin
	b_h = D_VTOL_prop + D_FW_prop
	S_h = horizontal_tail_area(V_H, mac_w, S_w, 0.75)
	c_h = rect_chord(S_h, b_h)
end 

# ╔═╡ 3519da90-3263-11eb-0f1b-392714a78701
begin
	x_PF = x_w + b_h / 2 * tan(Λ_LE_w) - (D_VTOL_prop / 2 - c_prop) / cos(Λ_LE_w)
	x_PR = x_w + c_r_w + b_h / 2 * tan(Λ_TE_w) + (D_VTOL_prop / 2 + c_prop) / cos(Λ_TE_w)
end

# ╔═╡ 40902c30-3268-11eb-32ef-75a5c60b85ab
md"""
Vertical tail.
"""

# ╔═╡ 8d97a350-30d3-11eb-307f-cd48a69b8300
V_V = 0.04

# ╔═╡ 03ad8000-3093-11eb-1068-19d569a57d29
begin
	x_CG = 0.5 * (x_PR - x_PF)
	x_VT = x_PR + D_VTOL_prop / 2 + c_prop
end

# ╔═╡ ff71bd40-3262-11eb-3155-9f4c31973837
begin	
	l_v = x_VT - x_CG
	λ_v = 0.5
	Λ_LE_v = deg2rad(20)
	S_v = vertical_tail_area(V_V, 1.0, S_w, l_v)
	b_v = vertical_tail_span(S_v, c_h, 0.8)
	AR_v = aspect_ratio(b_v, S_v)
	c_r_v = 2 * S_v / (1 + λ_v) / b_v
	c_t_v = λ_v * c_r_v
	mac_v = mean_aerodynamic_chord(c_r_v, λ_v)
	l_h = htail_arm(x_VT, b_v, Λ_LE_v, c_h, x_CG)
	l_v = vtail_arm(x_VT, mac_v)
end

# ╔═╡ 18c1fe50-3091-11eb-060c-c944f6e2e2c3
begin
    import Pkg
    Pkg.add(url="https://github.com/Pocket-titan/DarkMode.git")
    import DarkMode
	DarkMode.enable(theme="material-darker")
end

# ╔═╡ Cell order:
# ╟─23137140-3091-11eb-3f8a-2b8bd4b426b5
# ╟─cbac90c0-3262-11eb-1908-87d130f6d01c
# ╟─5fb5d550-3093-11eb-0691-09050e6939eb
# ╟─df735200-3263-11eb-3c15-c9f42ce9fbef
# ╠═304ebdb0-3091-11eb-0603-dd2131ef1226
# ╠═f5ae15b2-3091-11eb-3c08-53653717475a
# ╠═d5b5ec00-3092-11eb-2d10-87fcea8784a9
# ╟─8cfe52f0-3267-11eb-2b11-23458f7f2966
# ╠═7f48f152-3263-11eb-36d1-cfd2a24b2cd9
# ╟─98a56ace-3267-11eb-3627-ed82836734c2
# ╠═3d4d2870-3263-11eb-3e80-1f3f5db82243
# ╠═3519da90-3263-11eb-0f1b-392714a78701
# ╟─39d637e0-3268-11eb-2536-033d1b0d580e
# ╠═f11e9b60-3266-11eb-0c1f-fd0eaea6cc5d
# ╠═5c3009f0-3093-11eb-0ce2-9de2e3be4541
# ╟─40902c30-3268-11eb-32ef-75a5c60b85ab
# ╠═8d97a350-30d3-11eb-307f-cd48a69b8300
# ╠═03ad8000-3093-11eb-1068-19d569a57d29
# ╠═ff71bd40-3262-11eb-3155-9f4c31973837
# ╟─18c1fe50-3091-11eb-060c-c944f6e2e2c3
