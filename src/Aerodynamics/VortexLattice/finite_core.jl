# Finite-core model
function bound_leg_velocity(a, b, Γ, ε) 
    na, nb, σ = norm(a), norm(b), dot(a, b)
    term_1    = (na^2 - σ) / √(na^2 + ε^2) + (nb^2 - σ) / √(nb^2 + ε^2)
    term_2    = a × b / (na^2 * nb^2 - σ^2 + ε^2 * (na^2 + nb^2 - 2 * na * nb))
    
    Γ/4π * term_1 * term_2
end

trailing_leg_velocity(r, Γ, u, ε) = Γ/4π * normalize(r) × u / (norm(r) - dot(r, u) + ε^2 / (norm(r) + dot(r, u)))
trailing_legs_velocities(a, b, Γ, u, ε) = trailing_leg_velocity(a, Γ, u, ε) - trailing_leg_velocity(b, Γ, u, ε)
total_horseshoe_velocity(a, b, Γ, u, ε) = bound_leg_velocity(a, b, Γ, ε) + trailing_legs_velocities(a, b, Γ, u, ε)

horseshoe_velocity(r, line :: Line, Γ, direction, ε) = total_horseshoe_velocity(r1(r, line), r2(r, line), Γ, direction, ε)