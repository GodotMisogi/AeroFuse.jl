### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# ╔═╡ 55c656d0-3b7b-11eb-27c3-f52254894aec
md"""
# Constraints Analysis
"""

# ╔═╡ 993c5b40-406b-11eb-3d8b-356e309163f0
md"""
The constraints of the initial sizing procedure are determined by equating the thrust-to-weight ratio to a function of the wing loading for different flight conditions.
```math
(T/W) = f(W/S)
```
"""

# ╔═╡ d99d87f0-3b7e-11eb-30cb-e9db3f912b69
wing_loadings = 20:1:150;

# ╔═╡ 6ecda840-3b7b-11eb-2f36-a37e3aa55916
md"""
### Power Loading

\begin{align}
	(P/W) = \frac{(T/W) V}{\eta_\text{prop}}
\end{align}
"""

# ╔═╡ 3b52cb30-3b7b-11eb-1562-9f528047143b
power_loading(tbW, V, η_p) = tbW * V / η_p;

# ╔═╡ c25785d0-4071-11eb-23e1-b59dbe97ea29
η_prop_FW = 0.8;

# ╔═╡ 4ecf1920-3b80-11eb-063e-57d0c44af555
md"""
### Cruise
"""

# ╔═╡ 61cd2ad0-4071-11eb-197e-ebd10d747a99
md""" 
```math
\begin{align}
	(T/W)^\text{FW}_\text{cruise} = \frac{q C_{D_0}}{(W/S)} + \frac{k}{q}(W/S)
\end{align}
```
"""

# ╔═╡ 6beccf70-3b7b-11eb-12e6-d100322561d0
thrust_to_weight_fw_cruise(q, k, CD0, wing_loading) = q * CD0 / (wing_loading) + k / q * wing_loading;

# ╔═╡ 1b50c100-3b7d-11eb-1886-f7c9833720da
dynamic_pressure(ρ, V) = 1/2 * ρ * V^2;

# ╔═╡ 9e25b9d0-3b7f-11eb-1588-2d40cc56f40c
span_efficiency_factor(e, AR) = 1 / (π * e * AR);

# ╔═╡ 84a3c4b0-3b80-11eb-3819-09f88564eccc
begin
	cd0 = 0.004
	AR = 6
	e = 0.98
	V_cruise = 40
	ρ = 1.225
	k = span_efficiency_factor(e, AR)
	q = dynamic_pressure(ρ, V_cruise)
	
	tbWs_cruise = thrust_to_weight_fw_cruise.(q, k, cd0, wing_loadings)
	pbWs_cruise = power_loading.(tbWs_cruise, V_cruise, η_prop_FW)
end;

# ╔═╡ 411af330-3b80-11eb-17c6-611dfa0feb1e
md"""
### Rate of Climb
"""

# ╔═╡ f3ce6010-4072-11eb-170d-5b6f24a452d0
md"""
```math
V_\text{R/C} = \sqrt{\frac{2}{\rho} \left(\frac{W}{S}\right) \sqrt{\frac{k}{3 C_{D_0}}}}
```
"""

# ╔═╡ 8c993880-3b7b-11eb-07b4-51c590a54a38
best_climb_rate(wing_loading, k, CD0, ρ = 1.225) = sqrt((2 / ρ) * wing_loading * sqrt(k / 3CD0));

# ╔═╡ 0da45f10-3b84-11eb-2750-15deba04fc1e
V_RoC = best_climb_rate.(wing_loadings, k, cd0, ρ);

# ╔═╡ f0da3940-3b83-11eb-28fc-c144bad9a9a3
md"""
Fixed-wing
"""

# ╔═╡ 388948f0-4073-11eb-2396-6dd666582a97
md"""
```math
(T/W)^\text{FW}_\text{climb} = \frac{R/C}{V_\text{ROC}} + \frac{q C_{D_0}}{(W/S)} + \frac{k}{q}(W/S)
```
"""

# ╔═╡ 0b8b51b0-5570-11eb-0f19-7bf6626ebe1a
thrust_to_weight_fw_climb(wing_loading, rate_of_climb, best_climb_rate, q, CD0) = rate_of_climb / best_climb_rate + (q / wing_loading) * CD0 + (k / q) * wing_loading;

# ╔═╡ 65be11e0-3b80-11eb-1fab-0f96a6f46162
begin
	RoC = 3.

	tbWs_RoC = thrust_to_weight_fw_climb.(wing_loadings, RoC, V_RoC, q, cd0)
	pbWs_fw_climb = power_loading.(tbWs_RoC, V_RoC, η_prop_FW)
end;

# ╔═╡ fa9c6480-3b83-11eb-0dcc-7b93769f018f
md"""
VTOL
"""

# ╔═╡ 90406c90-4073-11eb-3217-7d0291296c4f
md"""
```math
(T/W)^\text{VTOL}_{climb} = 1.2\left(1 + \frac{1}{(W/S)}\rho (R/C)^2(S_\text{proj}/S_w)\right) 
```
"""

# ╔═╡ 91977310-3b7b-11eb-1391-67b38f3aaa31
thrust_to_weight_vtol_climb(wing_loading, V_vtol, S_proj_S_wing, CD, ρ = 1.225) = 1.2(1 + ρ * V_vtol^2 * S_proj_S_wing / wing_loading);

# ╔═╡ 41130720-3b84-11eb-2f38-fde35179f361
begin
	RoC_vtol = 2.
	CD = 2.
	S_proj_S_wing = 1.4
	
	tbWs_RoC_vtol = thrust_to_weight_vtol_climb.(wing_loadings, RoC_vtol, S_proj_S_wing, CD, ρ)
	pbWs_vtol_climb = power_loading.(tbWs_RoC_vtol, V_RoC, η_prop_FW)
end;

# ╔═╡ 44f60a80-3b80-11eb-08f1-c7e3d84c276d
md"""
### Stall Speed
"""

# ╔═╡ 058af970-4074-11eb-3b86-5766ff193840
md"""
```math
(W/S)^{FW}_\text{stall} = \frac{1}{2} \rho V_\text{stall} C_{L_\max}
```
"""

# ╔═╡ 25094f50-3b7d-11eb-2ab4-93c376e1b193
wing_loading_stall_speed(V_stall, CL_max, ρ = 1.225) = 1/2 * ρ * V_stall^2 * CL_max;

# ╔═╡ f6ac2640-3b7d-11eb-02de-21ef4c02ecac
begin
	V_stall = 10
	CL_max = 1.8
	
	n = 15
	stalls = fill(wing_loading_stall_speed(V_stall, CL_max), n);
	pbWs_stall = power_loading.(range(0, 2.; length = n), V_stall, CL_max)
end;

# ╔═╡ 5a681fe0-3b83-11eb-3da0-cbca7e455391
md"""
### Service Ceiling
"""

# ╔═╡ 621e7590-3b83-11eb-1700-91b2f71716d1
begin
	RoC_ceiling = 0.5
	
	tbWs_RoC_ceiling = thrust_to_weight_fw_climb.(wing_loadings, RoC_ceiling, V_RoC, q, cd0)
	pbWs_ceiling = power_loading.(tbWs_RoC_ceiling, V_RoC, η_prop_FW)
end;

# ╔═╡ 413bd390-3b7e-11eb-34f9-bf8ad1968b0f
begin
	using Plots
	plot(ylabel = "P/W", xlabel = "W/S", title = "Matching Chart (Fixed Wing)")
	plot!(stalls, pbWs_stall, label = "Stall Speed")
	plot!(wing_loadings, pbWs_cruise, label = "Cruise")
	plot!(wing_loadings, pbWs_fw_climb, label = "Fixed-Wing Climb")
	plot!(wing_loadings, pbWs_ceiling, label = "Ceiling")
end

# ╔═╡ 9a3c25c0-4070-11eb-25e2-05c4934b6332
md"""
## Matching Charts
"""

# ╔═╡ 824736f0-5572-11eb-2966-85bf07b88d2a
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ 67778190-5572-11eb-2db3-4d6202789a01
hint(md"""Call the `plot(x, y)` function""")

# ╔═╡ dd39bc4e-3b90-11eb-17df-09f9d206a6b2
begin
	plot(ylabel = "P/W", xlabel = "W/S", title = "Matching Chart (VTOL)")
	plot!(wing_loadings, pbWs_vtol_climb, label = "VTOL Climb")
end

# ╔═╡ Cell order:
# ╟─55c656d0-3b7b-11eb-27c3-f52254894aec
# ╟─993c5b40-406b-11eb-3d8b-356e309163f0
# ╠═d99d87f0-3b7e-11eb-30cb-e9db3f912b69
# ╟─6ecda840-3b7b-11eb-2f36-a37e3aa55916
# ╟─3b52cb30-3b7b-11eb-1562-9f528047143b
# ╟─c25785d0-4071-11eb-23e1-b59dbe97ea29
# ╟─4ecf1920-3b80-11eb-063e-57d0c44af555
# ╟─61cd2ad0-4071-11eb-197e-ebd10d747a99
# ╟─6beccf70-3b7b-11eb-12e6-d100322561d0
# ╟─1b50c100-3b7d-11eb-1886-f7c9833720da
# ╟─9e25b9d0-3b7f-11eb-1588-2d40cc56f40c
# ╟─84a3c4b0-3b80-11eb-3819-09f88564eccc
# ╟─411af330-3b80-11eb-17c6-611dfa0feb1e
# ╟─f3ce6010-4072-11eb-170d-5b6f24a452d0
# ╟─8c993880-3b7b-11eb-07b4-51c590a54a38
# ╟─0da45f10-3b84-11eb-2750-15deba04fc1e
# ╟─f0da3940-3b83-11eb-28fc-c144bad9a9a3
# ╟─388948f0-4073-11eb-2396-6dd666582a97
# ╟─0b8b51b0-5570-11eb-0f19-7bf6626ebe1a
# ╟─65be11e0-3b80-11eb-1fab-0f96a6f46162
# ╟─fa9c6480-3b83-11eb-0dcc-7b93769f018f
# ╟─90406c90-4073-11eb-3217-7d0291296c4f
# ╟─91977310-3b7b-11eb-1391-67b38f3aaa31
# ╟─41130720-3b84-11eb-2f38-fde35179f361
# ╟─44f60a80-3b80-11eb-08f1-c7e3d84c276d
# ╟─058af970-4074-11eb-3b86-5766ff193840
# ╟─25094f50-3b7d-11eb-2ab4-93c376e1b193
# ╟─f6ac2640-3b7d-11eb-02de-21ef4c02ecac
# ╟─5a681fe0-3b83-11eb-3da0-cbca7e455391
# ╟─621e7590-3b83-11eb-1700-91b2f71716d1
# ╟─9a3c25c0-4070-11eb-25e2-05c4934b6332
# ╠═824736f0-5572-11eb-2966-85bf07b88d2a
# ╟─67778190-5572-11eb-2db3-4d6202789a01
# ╟─413bd390-3b7e-11eb-34f9-bf8ad1968b0f
# ╟─dd39bc4e-3b90-11eb-17df-09f9d206a6b2
