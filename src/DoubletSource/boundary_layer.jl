module ViscFoil

include("DoubletSource.jl")
using .DoubletSource: doublet_matrix, source_matrix, source_strengths
using NLsolve

# function integral_momentum_equation(dx, x, p, s)
#     θ, dθ = x                       # State vector
#     ρ_e, u_e, a_e, τ_w, δ_st = p    # Parameters
#     H = δ_st/θ                      # Shape parameter
#     c_f = τ_w/(1/2 * ρ_e * u_e^2)   # Skin-friction coefficient
#     M_e, H = u_e/a_e                # Mach number
#     dx[1] = dθ                         # LHS: dθ/ds
#     dx[2] = c_f / 2  - (H + 2 - M_e^2) * θ/u_e * ??? # RHS: c_f/2  - (H + 2 - Mₑ²) * θ/U_e * dU_e/ds
    
# end

inviscid_edge_speed(panels :: AbstractVector{Panel2D}, uniform :: Uniform2D, source_influences, doublet_influences, ms, μs) = [ dot(velocity(uniform), panel_tangent(panel)) for panel in panels ] .+ source_influences * ms .+ doublet_influences * μs

δ_star(m, U_e) = m / U_e # δ* = m/U_e

shape_parameter(m, U_e, θ) = δ_star(m, U_e) / θ # H = δ*/θ

f1(H) = 1.515 + (H < 4 ? 0.076 : 0.040) * (H - 4)^2 / H
f2(H) = -0.067 + (H < 7.4 ? 0.01977 * (7.4 - H)^2 / (H - 1) : 0.022 * (1 - 1.4 / (H - 6))^2 )
f3(H) = 0.207 + (H < 4 ? 0.00205 * (4 - H)^5.5 : -0.003 * (H - 4)^2 )

cf(U_e, θ, H, ν = 1.5e-5) = 2ν / (U_e * θ) * f2(H)
c_Δ_times_2_by_H_star(ν, U_e, θ, H) = ν / (U_e * θ) * f3(H) # 2C_Δ/H* = f₃(H)/Re

R_1(θ, θ_1, U_e, U_e_1, x, x_1, H, c_f) = (θ - θ_1)/( (θ + θ_1) / 2) + (H + 2)*(U_e - U_e_1)/( (U_e + U_e_1) / 2) - (c_f / 2) * (x - x_1) / ( (θ + θ_1) / 2)

R_2(H_star, H_star_1, U_e, U_e_1, θ, θ_1, x, x_1, H, c_f, c_Δ) = (H_star - H_star_1)/( (H_star + H_star_1) / 2) + (1 - H) * (U_e - U_e_1) / ( (U_e + U_e_1) / 2) + (c_f / 2 - 2 * c_Δ / ((H_star + H_star_1) / 2)) * (x - x_1) / ((θ + θ_1) / 2)

Q(source_AIC, doublet_AIC, ms, μs, freestream_tangents) = source_AIC * ms .+ doublet_AIC * μs - freestream_tangents

function solve_system!(R, x, source_AIC, doublet_AIC, boco, )
    μ, Θ, m = x
    R1 = R_1()
    R2 = R_2()
    Q  = Q(source_AIC, doublet_AIC, ms, μs, boco)

    R .= [R1, R2, Q]
end