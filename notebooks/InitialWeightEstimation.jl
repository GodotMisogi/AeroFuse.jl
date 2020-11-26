### A Pluto.jl notebook ###
# v0.12.11

using Markdown
using InteractiveUtils

# ╔═╡ 77b76e80-2f3e-11eb-01d1-057a5351f181
using Pkg; Pkg.add("NLsolve")

# ╔═╡ b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
mass_takeoff(m_vtol_prop, m_fixed_prop, m_payload, ff_batt, ff_struct, ff_subsys, ff_avionics) = (m_vtol_prop + m_fixed_prop + m_payload) / (1 - (ff_batt + ff_struct + ff_subsys + ff_avionics))

# ╔═╡ Cell order:
# ╠═b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
# ╠═77b76e80-2f3e-11eb-01d1-057a5351f181
