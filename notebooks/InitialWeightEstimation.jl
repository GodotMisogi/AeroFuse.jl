### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
mass_takeoff(m_vtol_prop, m_fixed_prop, m_payload, ff_batt, ff_struct, ff_subsys, ff_avionics) = (m_vtol_prop + m_fixed_prop + m_payload) / (1 - (ff_batt + ff_struct + ff_subsys + ff_avionics))

# ╔═╡ 98de4950-3581-11eb-0495-05b7d0c5c57c
mass_residual(m_to, m_vtol_prop, m_fixed_prop, m_payload, mf_batt, mf_struct, mf_subsys, mf_avionics) = mto * (1 - (mf_batt + mf_struct + mf_subsys + mf_avionics)) - (m_vtol_prop + m_fixed_prop + m_payload)

# ╔═╡ 6e5a5bc0-3b75-11eb-27d1-dd09cf16244c
# mass_residual(40, _, _, 10, 0.3, 0.07)

# ╔═╡ 7bf733b0-3b76-11eb-1cfc-6f7225ac79d4
mf_battery(t, m_to, E_spec, η_batt, f_usable) = t * P / (m_to * E_spec * η_batt * f_usable)

# ╔═╡ 95f0ad50-3b76-11eb-03c7-ad71559db97b
mf_hover_battery(t_hover, E_spec, FM, η_elec, η_batt, f_usable, disk_loading, g = 9.81, ρ = 1.225) = t_hover * g / (E_spec * FM * η_elec * η_batt * f_usable) * √(disk_loading / 2ρ)

# ╔═╡ e406e950-3b76-11eb-1ea2-15bd92e22b1f


# ╔═╡ Cell order:
# ╠═b7cd7250-2f3c-11eb-260d-1f81b64e7cbe
# ╠═98de4950-3581-11eb-0495-05b7d0c5c57c
# ╠═6e5a5bc0-3b75-11eb-27d1-dd09cf16244c
# ╠═7bf733b0-3b76-11eb-1cfc-6f7225ac79d4
# ╠═95f0ad50-3b76-11eb-03c7-ad71559db97b
# ╠═e406e950-3b76-11eb-1ea2-15bd92e22b1f
