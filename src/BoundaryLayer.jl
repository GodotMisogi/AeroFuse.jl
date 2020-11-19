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

inviscid_edge_speed(panels :: Array{Panel2D}, uniform :: Uniform2D, source_influences, doublet_influences, ms, μs) = [ dot(velocity(uniform, panel_tangent(panel)) for panel in panels ] .+ source_influences * ms .+ doublet_influences * μs

δ_star(m, Uₑ) = m / Uₑ # δ* = m/Uₑ

shape_parameter(m, Uₑ, θ) = δ_star(m, Uₑ) / θ # H = δ*/θ

f1(H) = 1.515 + (H < 4 ? 0.076 : 0.040) * (H - 4)^2 / H
f2(H) = -0.067 + (H < 7.4 ? 0.01977 * (7.4 - H)^2 / (H - 1) : 0.022 * (1 - 1.4 / (H - 6))^2 )
f3(H) = 0.207 + (H < 4 ? 0.00205 * (4 - H)^5.5 : -0.003 * (H - 4)^2 )

cf(Uₑ, θ, H, ν = 1.5e-5) = 2ν / (Uₑ * θ) * f2(H)
c_Δ_times_2_by_H_star(ν, Uₑ, θ, H) = ν / (Uₑ * θ) * f3(H) # 2C_Δ/H* = f₃(H)/Re

R₁(θ, θ₁, Uₑ, Uₑ₁, x, x₁, H, c_f) = (θ - θ₁)/( (θ + θ₁) / 2) + (H + 2)*(Uₑ - Uₑ₁)/( (Uₑ + Uₑ₁) / 2) - (c_f / 2) * (x - x₁) / ( (θ + θ₁) / 2)
R₂(H_star, H_star₁, Uₑ, Uₑ₁, θ, θ₁, x, x₁, H, c_f, c_Δ) = (H_star - H_star₁)/( (H_star + H_star₁) / 2) + (1 - H) * (Uₑ - Uₑ₁) / ( (Uₑ + Uₑ₁) / 2) + (c_f / 2 - 2 * c_Δ / ((H_star + H_star₁) / 2)) * (x - x₁) / ((θ + θ₁) / 2)
R₃(source_influences, doublet_influences, ms, μs, freestream_tangents) = source_influences * ms .+ doublet_influences * μs - freestream_tangents 

end