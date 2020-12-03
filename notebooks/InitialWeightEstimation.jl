### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
mass_takeoff(m_vtol_prop, m_fixed_prop, m_payload, ff_batt, ff_struct, ff_subsys, ff_avionics) = (m_vtol_prop + m_fixed_prop + m_payload) / (1 - (ff_batt + ff_struct + ff_subsys + ff_avionics))

# ╔═╡ 98de4950-3581-11eb-0495-05b7d0c5c57c


# ╔═╡ 3f5de762-3260-11eb-071d-6d2fb6839e0e
begin
    # import DarkMode
	# DarkMode.enable(theme="material-darker")
end

# ╔═╡ Cell order:
# ╠═b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
# ╠═98de4950-3581-11eb-0495-05b7d0c5c57c
# ╟─3f5de762-3260-11eb-071d-6d2fb6839e0e
