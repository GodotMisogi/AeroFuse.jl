# Generic
dynamic_pressure(ρ, V) = 0.5 * ρ * V^2
stall_speed(wing_loading, CL_max, ρ = 1.225) = 2 * wing_loading/(ρ * CL_max)

# Drag polar
span_efficiency_factor(e, AR) = 1 / (π * e * AR)
drag_polar(CD0, k, CL) = CD0 + k * CL^2

# Wing loading
wing_loading_stall_speed(V_stall, CL_max, ρ = 1.225) = 0.5 * ρ * V_stall^2 * CL_max


## Weight estimation and performance
#==========================================================================================#

# Power loading
power_loading(tbW, V, η_p) = tbW * V / η_p

# Thrust-to-weight ratios
thrust_to_weight_fw_cruise(q, CD0, wing_loading) = q * CD0 * 1 / (wing_loading) + k / q * wing_loading

best_climb_rate(wing_loading, k, CD0, ρ = 1.225) = sqrt((2 / ρ) * wing_loading * sqrt(k / (3 * CD0))) 
thrust_to_weight_fw_climb(rate_of_climb, best_climb_rate, q, CD0) = rate_of_climb / best_climb_rate + q / wing_loading * CD0 + k / q * wing_loading

thrust_to_weight_vtol_climb(wing_loading, V_vtol, S_total, S_wing, CD, ρ = 1.225) = 1.2(1 + ρ * V_vtol^2 * S_total / S_wing / wing_loading)
# W + 0.5 * ρ * V_vtol^2 * S_proj * CD

# Mass estimation
mass_takeoff(m_vtol_prop, m_fixed_prop, m_payload, ff_batt, ff_struct, ff_subsys, ff_avionics) = (m_vtol_prop + m_fixed_prop + m_payload) / (1 - (ff_batt + ff_struct + ff_subsys + ff_avionics))

flat_plate_cd(α) = 2(sin(α))^3

figure_of_merit(thrust) = 0.4742 * thrust^0.0793

hover_speed(T, S_p, ρ) = sqrt(0.5 * T / (ρ * S_p))

induced_velocity(climb_rate_vtol, v_h) = let x = climb_rate_vtol / (2v_h); -x + sqrt(x^2 + 1) end

disk_loading(M_TO) = 3.2261M_TO + 74.991

prop_area(W_TO, DL, η_prop) = W_TO / (DL * η_prop)

area(a, b, c, d) = a * b * c / d
prodsratio(a, b, c, d) = a * b / (c * d) 