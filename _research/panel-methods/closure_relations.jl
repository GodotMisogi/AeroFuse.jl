## CLOSURE RELATIONS
#============================================#

# Momentum shape for density: H
H_Whitfield(H, M_e) = (H - 0.290M_e^2)/(1 + 0.113M_e^2)

## Laminar closure for walls - Falkner-Skan profiles
#============================================#

# Laminar kin-ergy shape: H*
lam_H_star(H_k) = 1.515 + ifelse(H_k < 4., 0.076, 0.040) * (H_k - 4)^2 / H_k

# Laminar skin-friction: Re_θ CF/2
lam_CF(H_k) = -0.067 + ifelse(H_k < 7.4, 0.01977(7.4 - H_k)^2 / (H_k - 1), 0.022(1. - 1.4 / (H_k - 6.))^2)

# Laminar dissipation: Re_θ 2CD / H_k*
lam_CD(H_k) = 0.207 + ifelse(H_k < 4., 0.00205(4 - H_k)^5.5, -0.003(H_k - 4.)^2 / (1. + 0.02(H_k - 4.)^2))

## Laminar closure for wakes
#============================================#

# Laminar wake kin-ergy shape: H*
lam_wake_H_star(H_k) = 1.50 + ifelse(H_k < 3.5, 0.025(3.5 - H_k)^3 + 0.001(3.5 - H_k)^5, 0.015(H_k - 3.5)^2 / H_k)

# Laminar wake dissipation: Re_θ 2CD / H*
lam_wake_CD(H_k) = 1.52(H_k - 1)^2 / (3 + H_k^3)


## Turbulent closure relations
#============================================#

# Skin-friction
turb_CF(H_k, Re_θ, M_e) = let F_c = √(1 + 0.2M_e^2); 0.3exp(-1.33H_k) * log(10, Re_θ / F_c)^(-1.74 - 0.31H_k) / F_c + 1.1e-4(tanh(4 - H_k/0.875) - 1) end

# Density thickness (negligible at low speeds): H**
turb_H_star_star(H_k, M_e) = (0.064 / (H_k - 0.8) + 0.251)M_e^2

# Y+
y_plus(u_t, μ_e, ρ_e, η) = ρ_e * u_t * η / μ_e

# Velocity profile
function turb_ue(Cf, θ, y_p, η, a, b)
	uτ_ue = √(Cf / 2)
	s 	  = Cf / abs(Cf)

	uτ_ue * s / 0.09 * atan(0.09y_p) + (1 - uτ_ue * s * π / 0.18) * √(tanh(a * (η, θ)^b))
end

# H*?
function turb_H_star(H_k, Re_θ)
	H_0 = ifelse(Re_θ < 400, 4, 3 + 400/Re_θ)

	1.505 + 4/Re_θ + ifelse(H_k < H_0, (0.165 - 1.6/√Re_θ) * (H_0 - H_k)^1.6 / H_k, (H_k - H_0)^2 * (0.04 / H_k + 0.007log(10, Re_θ)/(H_k - H_0 + 4/log(10, Re_θ))^2))
end

# CD
turb_CD(cf, u_s, c_τ) = cf * u_s / 2 + c_τ * (1 - u_s)

# U_s
wall_slip(H_star, H_k, H) = H_star / 2 * (1 - 4/3(H_k - 1) / H)

# Equilibirum shear stress coefficient C_τ_EQ
c_τ_EQ(H, H_k, H_star, u_s) = 0.015H_star / (1. - u_s) * (H_k - 1.)^3 / (H_k^2 * H)


## shape parameters for separation criteria
hlmax = 3.8;
htmax = 2.5;


## Transition
#============================================#

# Rate slope dñ/dReθ
dñdReθ(H_k) = 0.01 * √((2.4H_k - 3.7 + 2.5tanh(1.5H_k - 4.65))^2 + 0.25)

# Critical Reynolds number: log₁₀ Re_θ
log10Reθ(H_k) = (1.415/(H_k - 1) - 0.489)tanh(20/(H_k - 1) - 12.9) + 3.295/(H_k - 1) + 0.44

# G-β locus?
G(β) = 6.7 * √(1 + 0.75β)

# Weird amplification ratio DE
l(H_k) = (6.54H_k - 14.07) / H_k^2
m(H_k) = (0.058(H_k - 4)^2 / (H_k - 1) - 0.068) / l(H_k)



## Generic conversions
#============================================#

# Cf = 2lam_CF(H) / Re
cf(U_e, θ, H, ν = 1.5e-5) = 2 * lam_CF(H) / reynolds_number(U_e, θ, ν) 

# CD = lam_H_star(H) * lam_CD(H) / 2Re
cD(U_e, θ, H, ν) = lam_H_star(H) * lam_CD(H) / (2 * reynolds_number(U_e, θ, ν) )

