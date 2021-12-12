module ViscFoil

include("VorticityStream.jl")
# function integral_momentum_equation(dx, x, p, s)
#     θ, dθ = x                       # State vector
#     ρ_e, u_e, a_e, τ_w, δ_st = p    # Parameters
#     H = δ_st/θ                      # Shape parameter
#     c_f = τ_w/(1/2 * ρ_e * u_e^2)   # Skin-friction coefficient
#     M_e, H = u_e/a_e                # Mach number
#     dx[1] = dθ                      # LHS: dθ/ds
#     dx[2] = c_f / 2  - (H + 2 - M_e^2) * θ/u_e * ??? # RHS: c_f/2  - (H + 2 - Mₑ²) * θ/U_e * dU_e/ds
# end

function edge_velocity(σs, panels, wake_panels, u)
    all_panels          = [panels; wake_panels]
    foil_vortices       = vortex_matrix(panels, panels)
    all_sources         = source_matrix(panels, all_panels)
    wake_wake_vortices  = vortex_matrix(wake_panels, wake_panels)
    wake_foil_vortices  = vortex_matrix(wake_panels, panels)
    wake_boundary       = boundary_vector(colpoints.(wake_panels), u)

    P       = -foil_vortices^(-1) * all_sources
    P_s     = 
    # D1  = visc_visc_vortices # PLACEHOLDER
    D1      = wake_foil_vortices * P + wake_wake_vortices
    U_inv   = [ (zeros ∘ length)(panels) ; 
                    wake_boundary  ]
    D       = [ P
                D1 ]

    U_e     = U_inv + D * σs
end

end