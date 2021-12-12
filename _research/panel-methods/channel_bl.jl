using AeroMDAO

AeroMDAO.reynolds_number(U, L, ν) 	= U * L / ν

## Correlation functions

# Laminar kin-ergy shape: H*
lam_H_star(H) = 1.515 + ifelse(H < 4., 0.076, 0.040) * (H - 4)^2 / H

# Laminar skin-friction: Re_θ CF/2
lam_CF(H) = -0.067 + ifelse(H < 7.4, 0.01977(7.4 - H)^2 / (H - 1), 0.022(1. - 1.4 / (H - 6.))^2)

# Laminar dissipation: Re_θ 2CD / H*
lam_CD(H) = 0.207 + ifelse(H < 4., 0.00205(4 - H)^5.5, -0.003(H - 4.)^2 / (1. + 0.02(H - 4.)^2))

## Channel flow setup
u_edge(ρ_e, m_chan, h, δ_star) = m_chan / (ρ_e * (h - δ_star))

function channel_flow(h, H, δ_star, H_star, dHstar_dH)
    β_Hs_H = H / H_star * dHstar_dH
    B_δs_u = -δ_star / (h - δ_star)

    [   1       0   H + 2;
     -β_Hs_H β_Hs_H 1 - H;
        0    B_δs_u   1  ]
end

channel_sources(x, θ, δ_star, H_star, c_f, c_D, h, dh_dx) = [ x * c_f / 2θ; x / θ * (2c_D / H_star - c_f / 2); -x / (h - δ_star) * dh_dx ]

update_variables(p, β_p, x, x_prev) = p * (x / x_prev)^β_p

function solve_channel(ρ :: Real, h1 :: Real, u_invs :: Vector{<: Real}, xs :: Vector{<: Real})
    L = sum(xs)
    ν = 1.5e-5

    m_chan = ρ * u_invs[1] * h1

    hs = h1 * u_invs[1] ./ u_invs

    θ, δ = 0.2, 0.3

    solve_channel(hs, xs, u_invs, θ, δ, ν)
end

function solve_channel(m_chan :: Real, hs :: Vector{<: Real}, xs :: Vector{<: Real}, ρ :: Real)
    L = sum(xs)
    ν = 1.5e-5

    u_invs = m_chan / ρ ./ hs

    θ, δ = 0.1, 0.2

    solve_channel(hs, xs, u_invs, θ, δ, ν :: Real)
end

function compute_params(δ, θ, u_e, ν)
    H       = δ / θ
    Re_θ    = reynolds_number(u_e, θ, ν)
    H_star  = lam_H_star(H)
    c_f     = lam_CF(H) * 2 / Re_θ
    c_D     = lam_CD(H) * 2H_star / Re_θ
    H, Re_θ, H_star, c_f, c_D
end

function solve_system(h, H, δ_star, H_star, dHstar_H, x, θ, c_f, c_D, dhdx)
    A = channel_flow(h, H, δ_star, H_star, dHstar_H)
    b = channel_sources(x, θ, δ_star, H_star, c_f, c_D, h, dhdx)

    βs = A \ b
end

function solve_channel(hs :: Vector{<: Real}, xs :: Vector{<: Real}, u_invs :: Vector{<: Real}, θ :: Real, δ :: Real, ν :: Real)
    n = length(u_invs)
    u1 = first(u_invs)

    # State variables
    θs      = [ θ; zeros(n-1) ]
    δ_stars = [ δ; zeros(n-1) ]
    u_es    = [ u1; zeros(n-1) ]

    # State parameters
    H, Re_θ, H_star, c_f, c_D = compute_params(δ, θ, u1, ν)
    Hs      = [ H; zeros(n-1) ]
    H_stars = [ H_star; zeros(n-1) ]
    Re_θs   = [ Re_θ; zeros(n-1) ]
    c_fs    = [ c_f; zeros(n-1) ]
    c_Ds    = [ c_D; zeros(n-1) ]

    dHstar_H    = H / H_star
    dhdx        = 0.
    for i in 2:n
        βs          = solve_system(hs[i-1], Hs[i-1], δ_stars[i-1], H_stars[i-1], dHstar_H, xs[i-1], θs[i-1], c_fs[i-1], c_Ds[i-1], dhdx)

        θs[i]       = update_variables(θs[i-1], βs[1], xs[i], xs[i-1])
        δ_stars[i]  = update_variables(δ_stars[i-1], βs[2], xs[i], xs[i-1])
        u_es[i]     = update_variables(u_es[i-1], βs[3], xs[i], xs[i-1])

        Hs[i], Re_θs[i], H_stars[i], c_fs[i], c_Ds[i] = compute_params(δ, θ, u1, ν)
        dHstar_H    = (H_stars[i] - H_stars[i-1]) / (Hs[i] - Hs[i-1])
        dhdx        = (hs[i] - hs[i-1]) / (xs[i] - xs[i-1])
    end

    hs, θs, δ_stars, u_invs, u_es, Hs, H_stars, c_Ds, c_fs, Re_θs
end

## Case 1: U_invs specified
xs = collect(range(0, 1.0, length = 50))
u_invs = exp.(-xs)
h = 0.5
hs, θs, δ_stars, u_invs, u_es, Hs, H_stars, c_Ds, c_fs, Re_θs = solve_channel(1.225, h, u_invs, xs)

##
using Plots
plot(xs, u_invs, label = "U_invs")
plot!(xs[1:end-1], u_es, label = "U_es")
plot!(xs, hs ./ h, label = "h/h_inlet")

##
plot(xs[1:end-1], Hs[1:end-1])

## Case 2: hs specified
xs = collect(range(0, 1.0, length = 50))
h0, h1 = 0.1, 4.
hs = h0 .+ (h1 - h0) .* ((1.35 .* xs) ./ (0.35 .+ xs)).^2
m_chan = 0.1
hs, θs, δ_stars, u_invs, u_es, Hs, H_stars, c_Ds, c_fs, Re_θs = solve_channel(m_chan, hs, xs, 1.225)

##
# plot(xs[1:end-1], u_invs, label = "U_invs")
# plot!(xs[1:end-1], u_es, label = "U_es")
plot(xs, hs, label = "h")

##
plot(xs[1:end-1], Hs[1:end-1])
