"""
Reference: Drela's QPROP - https://web.mit.edu/drela/Public/web/qprop/qprop_theory.pdf
"""

function qprop_airfoil_lift_coefficient_model(α, Ma)
    alpha_rad = deg2rad(α)
    beta      = (1 - Ma^2)^0.5
    cl_0      = 0.5
    cl_alpha  = 5.8
    cl_min    = -0.3
    cl_max    = 1.2
    cl        = (alpha_rad * cl_alpha + cl_0) / beta
    Cl        = minimum(maximum(cl, cl_min), cl_max)
end

averaged_tangential_velocity(F, λ_w, B, r) = vt * F * √(1 + (4λ_w * R / (π * B * r))^2)
local_tangential_velocity(B, Γ, F, r, R, λ_w) = B * Γ / (4π * r) * 1 / (F * √(1 + (4λ_w * R / (π * B * r))^2))
local_wake_advance_ratio(r, R, φ) = r / R * tan(φ)
blade_circulation(W, c, cl) = 0.5 * W * c * cl

function prandtl_factor(B, r, R, φ)
    λ_w = local_wake_advance_ratio(r, R, φ)
    f = B / 2 * (1 - r / R) * 1 / λ_w
    F = 2 / π * (acos ∘ exp)(-f)
end