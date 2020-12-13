### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 55c656d0-3b7b-11eb-27c3-f52254894aec
md"""
# Constraints Analysis
"""

# ╔═╡ 1b50c100-3b7d-11eb-1886-f7c9833720da
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2

# ╔═╡ 9e25b9d0-3b7f-11eb-1588-2d40cc56f40c
span_efficiency_factor(e, AR) = 1 / (π * e * AR)

# ╔═╡ 6ecda840-3b7b-11eb-2f36-a37e3aa55916
md"""
### Power Loading

\begin{align}
	(P/W) = \frac{(T/W) V}{\eta_\text{prop}}
\end{align}
"""

# ╔═╡ 3b52cb30-3b7b-11eb-1562-9f528047143b
power_loading(tbW, V, η_p) = tbW * V / η_p;

# ╔═╡ d99d87f0-3b7e-11eb-30cb-e9db3f912b69
begin
	cd0 = 0.004
	AR = 6
	e = 0.98
	V_cruise = 40
	η_prop_FW = 0.8
	wing_loadings = 20:1:150
	ρ = 1.225
	k = span_efficiency_factor(e, AR)
	q = dynamic_pressure(ρ, V_cruise)
	CL_max = 1.8
	V_stall = 10
end;

# ╔═╡ 4ecf1920-3b80-11eb-063e-57d0c44af555
md"""
### Cruise
"""

# ╔═╡ 6beccf70-3b7b-11eb-12e6-d100322561d0
thrust_to_weight_fw_cruise(q, k, CD0, wing_loading) = q * CD0 * 1 / (wing_loading) + k / q * wing_loading;

# ╔═╡ 84a3c4b0-3b80-11eb-3819-09f88564eccc
begin
	tbWs_cruise = thrust_to_weight_fw_cruise.(q, k, cd0, wing_loadings)
	pbWs_cruise = power_loading.(tbWs_cruise, V_cruise, η_prop_FW)
end;

# ╔═╡ 411af330-3b80-11eb-17c6-611dfa0feb1e
md"""
### Rate of Climb
"""

# ╔═╡ 8c993880-3b7b-11eb-07b4-51c590a54a38
best_climb_rate(wing_loading, k, CD0, ρ = 1.225) = sqrt((2 / ρ) * wing_loading * sqrt(k / (3 * CD0)));

# ╔═╡ 0da45f10-3b84-11eb-2750-15deba04fc1e
V_RoC = best_climb_rate.(wing_loadings, k, cd0, ρ);

# ╔═╡ f0da3940-3b83-11eb-28fc-c144bad9a9a3
md"""
Fixed-wing
"""

# ╔═╡ 919724f0-3b7b-11eb-35ce-bbb61bfb605f
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

# ╔═╡ 91977310-3b7b-11eb-1391-67b38f3aaa31
thrust_to_weight_vtol_climb(wing_loading, V_vtol, S_proj_S_wing, CD, ρ = 1.225) = 1.2(1 + ρ * V_vtol^2 * S_proj_S_wing / wing_loading);

# ╔═╡ 41130720-3b84-11eb-2f38-fde35179f361
begin
	RoC_vtol = 2.
	CD = 2.
	S_proj_S_wing = 1.2
	
	tbWs_RoC_vtol = thrust_to_weight_vtol_climb.(wing_loadings, RoC_vtol, S_proj_S_wing, CD, ρ)
	pbWs_vtol_climb = power_loading.(tbWs_RoC_vtol, V_RoC, η_prop_FW)
end;

# ╔═╡ 44f60a80-3b80-11eb-08f1-c7e3d84c276d
md"""
### Stall Speed
"""

# ╔═╡ 25094f50-3b7d-11eb-2ab4-93c376e1b193
wing_loading_stall_speed(V_stall, CL_max, ρ = 1.225) = 0.5 * ρ * V_stall^2 * CL_max;

# ╔═╡ f6ac2640-3b7d-11eb-02de-21ef4c02ecac
begin
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

# ╔═╡ dd39bc4e-3b90-11eb-17df-09f9d206a6b2
begin
	plot(ylabel = "P/W", xlabel = "W/S", title = "Matching Chart (VTOL)")
	plot!(wing_loadings, pbWs_vtol_climb, label = "VTOL Climb")
end

# ╔═╡ Cell order:
# ╟─55c656d0-3b7b-11eb-27c3-f52254894aec
# ╠═1b50c100-3b7d-11eb-1886-f7c9833720da
# ╠═9e25b9d0-3b7f-11eb-1588-2d40cc56f40c
# ╟─6ecda840-3b7b-11eb-2f36-a37e3aa55916
# ╠═3b52cb30-3b7b-11eb-1562-9f528047143b
# ╠═d99d87f0-3b7e-11eb-30cb-e9db3f912b69
# ╟─4ecf1920-3b80-11eb-063e-57d0c44af555
# ╠═6beccf70-3b7b-11eb-12e6-d100322561d0
# ╠═84a3c4b0-3b80-11eb-3819-09f88564eccc
# ╟─411af330-3b80-11eb-17c6-611dfa0feb1e
# ╠═8c993880-3b7b-11eb-07b4-51c590a54a38
# ╠═0da45f10-3b84-11eb-2750-15deba04fc1e
# ╟─f0da3940-3b83-11eb-28fc-c144bad9a9a3
# ╠═919724f0-3b7b-11eb-35ce-bbb61bfb605f
# ╠═65be11e0-3b80-11eb-1fab-0f96a6f46162
# ╟─fa9c6480-3b83-11eb-0dcc-7b93769f018f
# ╠═91977310-3b7b-11eb-1391-67b38f3aaa31
# ╠═41130720-3b84-11eb-2f38-fde35179f361
# ╟─44f60a80-3b80-11eb-08f1-c7e3d84c276d
# ╠═25094f50-3b7d-11eb-2ab4-93c376e1b193
# ╠═f6ac2640-3b7d-11eb-02de-21ef4c02ecac
# ╟─5a681fe0-3b83-11eb-3da0-cbca7e455391
# ╠═621e7590-3b83-11eb-1700-91b2f71716d1
# ╠═413bd390-3b7e-11eb-34f9-bf8ad1968b0f
# ╠═dd39bc4e-3b90-11eb-17df-09f9d206a6b2
