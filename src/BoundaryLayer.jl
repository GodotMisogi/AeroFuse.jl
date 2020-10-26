module BoundaryLayer

# function integral_momentum_equation(dx, x, p, s)
#     θ, dθ = x                       # State vector
#     ρ_e, u_e, a_e, τ_w, δ_st = p    # Parameters
#     H = δ_st/θ                      # Shape parameter
#     c_f = τ_w/(0.5 * ρ_e * u_e^2)   # Skin-friction coefficient
#     M_e, H = u_e/a_e                # Mach number
#     dx[1] = dθ                         # LHS: dθ/ds
#     dx[2] = c_f / 2  - (H + 2 - M_e^2) * θ/u_e * ??? # RHS: c_f/2  - (H + 2 - Mₑ²) * θ/uₑ * duₑ/ds
    
# end

inviscid_edge_speed(panels :: Array{Panel2D}, uniform :: Uniform2D, ms, μs) = [ dot(velocity(uniform, panel_tangent(panel)) for panel in panels ] .+ ms .+ μs

displacement_thickness(m, Uₑ) = m / Uₑ # δ* = m/Uₑ

shape_parameter(m, Uₑ, θ) = displacement_thickness(m, Uₑ) / θ # H = δ*/θ

# kinetic_energy_shape_parameter()

f1(H) = 1.515 + (H < 4 ? 0.076 : 0.040) * (H - 4)^2 / H
f2(H) = -0.067 + (H < 7.4 ? 0.01977 * (7.4 - H)^2 / (H - 1) : 0.022 * (1 - 1.4 / (H - 6))^2 )
f3(H) = 0.207 + (H < 4 ? 0.00205 * (4 - H)^5.5 : -0.003 * (H - 4)^2 )

cf_by_2(ν, Uₑ, θ, H) = ν / (Uₑ * θ) * f2(H)
c_Δ_times_2_by_H_star(ν, Uₑ, θ, H) = ν / (Uₑ * θ) * f3(H) 



end