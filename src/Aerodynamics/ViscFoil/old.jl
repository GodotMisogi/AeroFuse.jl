
	# props			= laminar_body_quantities.(H_ks, M_es, Re_θs)
	# H_stars 		= getindex.(props, 1)
	# H_star_stars 	= getindex.(props, 2)
	# CFs 			= getindex.(props, 3)
	# CDs 			= getindex.(props, 4)

	# # Forward differencing boundary layer inputs using logarithmic derivatives
	
	# # Logarithmic differences
	# dlogUs	 		= @. log(Us[2:end] / Us[1:end-1])
	# dlogθs 	 		= @. log(θs[2:end] / θs[1:end-1])
	# dlogH_stars 	= @. log(H_stars[2:end] / H_stars[1:end-1])
	
	# # Differences
	# ΔTs 			= @. Ts[2:end] - Ts[1:end-1]

	# # Points
	# M_e_sq_as 		= weighted_vector(M_es[2:end], M_es[1:end-1], 0.5).^2 # Check Mₑ² term later
	# θ_as 	 		= weighted_vector(θs[2:end], θs[1:end-1], 0.5)
	# H_as 	 		= weighted_vector(Hs[2:end], Hs[1:end-1], 0.5)
	# H_k_as 	 		= weighted_vector(H_ks[2:end], H_ks[1:end-1], 0.5)
	# CF_as 	 		= weighted_vector(CFs[2:end], CFs[1:end-1], 0.5)
	# CD_as 	 		= weighted_vector(CFs[2:end], CFs[1:end-1], 0.5)
	# H_star_as 	 	= weighted_vector(H_stars[2:end], H_stars[1:end-1], 0.5)
	# H_star_star_as 	= weighted_vector(H_star_stars[2:end], H_star_stars[1:end-1], 0.5)
	# δ_star_as 		= weighted_vector(δ_stars[2:end], δ_stars[1:end-1], 0.5)

	# # Momentum equation residuals
	# R1 		 		= momentum_fd.(dlogθs, dlogUs, ΔXs, θ_as, H_k_as, CF_as, M_e_sq_as)
	# # Shape parameter equation residuals
	# R2 		 		= shape_fd.(dlogH_stars, dlogUs, ΔXs, θ_as, H_k_as, H_star_as, H_star_star_as, CF_as, CD_as)
	# # eⁿ amplification residuals
	# R3 				= lam_turb.(Ts[1:end-1], Ts[2:end], ΔXs, dlogUs, θ_as, δ_star_as, H_as, H_k_as, H_star_as, CF_as, n_crit)
	
	# n 		= length(θs)
	# H_stars = zeros(n)
	# H_star_stars = zeros(n)
	# CFs 	= zeros(n)
	# CDs 	= zeros(n)
	# U_s 	= zeros(n)
	# c_τs 	= zeros(n)
	# for i in 1:n
	# 	if tags[i] <: Panel2D
	# 		# Need to put transition checks
	# 		H_stars[i], H_star_stars[i], CFs[i], CDs[i] = laminar_body_quantities(H_ks[i], M_es[i], Re_θs[i])
	# 	else # Wake panel condition
	# 		H_stars[i], H_star_stars[i], CFs[i], c_τs[i], U_s[i], CDs[i] = turbulent_body_quantities(H_ks[i], Hs[i], M_es[i], Re_θs[i])
	# 	end
	# end

	# Laminar/transition/turbulence residuals
	# dlogc_τs 		= @. log(c_τs[2:end] / c_τs[1:end-1])
	# c_τ_EQ_as 		= weighted_vector(c_τs[2:end], c_τs[1:end-1], 0.5)

	# for i in 1:n
	# 	if tags[i] <: Panel2D
	# 		# Need to put transition checks
	# 		R3[i] = amplification_fd(Δn, dlogs, H_as, θ_as)
	# 	else # Wake panel condition
	# 		R3[i] = stress_transport_fd.(c_τ_as, dlogc_τs, dlogUs, ΔXs, δ, δ_stars, CF_as, c_τ_EQ_as, H_as)
	# 	end
	# end