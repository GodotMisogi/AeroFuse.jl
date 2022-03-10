## Upwinding
#============================================#

function upwind_factor(Hk1, Hk2, wake = false)
    Hut = 1.0
    C   = ifelse(wake, 1.0, 5.0)
    Huc = C * Hut / Hk2^2
    aa  = (Hk2 - 1) / (Hk1 - 1)
    la  = log(sign(aa) * aa)
    Hls = ifelse(la^2 > 15., 15., la^2)
    upw = 1. - exp(-Hls * Huc) / 2

    return upw
end

upwind(upw, f1, f2) = (1 - upw) * f1 + upw * f2


## BOUNDARY LAYER SYSTEM (UNUSED)
#============================================#

struct BoundaryLayer2D{T <: Real}
    θ       :: T
    δs      :: T
    tu      :: T
    Ue      :: T
    H       :: T
    Hk      :: T
    Hs      :: T
    Hss     :: T
    cf      :: T
    cDi     :: T
    cτ_EQ   :: T
    Us      :: T
    Me      :: T
    Reθ     :: T
end

function BoundaryLayer2D(θ, δs, p, Ue, params; turb = false, wake = false)
    # Type promotion for automatic differentiation
    T = promote_type(eltype(θ), eltype(δs), eltype(p), eltype(Ue))

    # Unpack parameters
    γ, μ0, ρ0, h0, Tsrat, KT_β, KT_λ, Re, M∞, V∞, K_lag, GA, GB, GC, ηD, cτ_C, cτ_E, ñ_crit = params

    H   = δs / θ
    Uk  = 0.
    Me  = 0.
    Reθ = 0.
    Hk  = 0.
    Hss = 0.

    if M∞ > 0
        Uk  = karman_tsien_correction(Ue, V∞, KT_λ)
        Me  = mach_number(Uk, h0, γ)

        Reθ = reynolds_number(θ, Uk, M∞, h0, γ, Tsrat, μ0, ρ0)

        Hk  = whitfield_shape_parameter(H, Me)
        Hss = turbulent_density_shape_parameter(Hk, Me)
    else 
        Uk  = Ue
        Reθ = reynolds_number(ρ0, Uk, θ, μ0)
        Hk  = H
    end

    if turb 
        sqrt_cτ = p
        Hs, cf, cDi, cτ_EQ, Us = turbulent_body_quantities(sqrt_cτ^2, Hk, H, Me, Reθ, GA, GB, GC, wake)

        return BoundaryLayer2D(θ, δs, p, Uk, H, Hk, Hs, Hss, cf, cDi * Hs / 2, cτ_EQ, Us, Me, Reθ)
    else 
        Hs, cf, cDi = laminar_body_quantities(Hk, Reθ)

        cτ_EQ = 0.
        Us    = 0.

        return BoundaryLayer2D{T}(θ, δs, p, Uk, H, Hk, Hs, Hss, cf, cDi * Hs / 2, cτ_EQ, Us, Me, Reθ)
    end
end

BoundaryLayer2D(U, params; turb = false, wake = false) = BoundaryLayer2D(U[1], U[2], U[3], U[4], params; turb = turb, wake = wake)

# State vector
state(bl :: BoundaryLayer2D) = [bl.θ; bl.δs; bl.tu; bl.Ue]

# Upwinding stations
function upwind(bl1 :: BoundaryLayer2D, bl2 :: BoundaryLayer2D; wake = false)
    upw = upwind_factor(bl1.Hk, bl2.Hk, wake)
    BoundaryLayer2D(map(x -> upwind(upw, getfield(bl1, x), getfield(bl2, x)), fieldnames(typeof(bl1)))...)
end