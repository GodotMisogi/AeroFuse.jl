module BoundaryLayer

function integral_momentum_equation(dx, x, p, s)
    θ, dθ = x                       # State vector
    ρ_e, u_e, a_e, τ_w, δ_st = p    # Parameters
    H = δ_st/θ                      # Shape parameter
    c_f = τ_w/(0.5 * ρ_e * u_e^2)   # Skin-friction coefficient
    M_e, H = u_e/a_e                # Mach number
    dx[1] = dθ                         # LHS: dθ/ds
    dx[2] = c_f / 2  - (H + 2 - M_e^2) * θ/u_e * ??? # RHS: c_f/2  - (H + 2 - Mₑ²) * θ/uₑ * duₑ/ds
    
end

end