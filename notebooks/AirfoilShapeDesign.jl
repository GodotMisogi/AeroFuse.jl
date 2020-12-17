### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# ╔═╡ 9c8e3e2e-4093-11eb-3e43-85d1fead00e5
begin
	using Revise
	using AeroMDAO
	using Optim
	using Plots
	gr()
end

# ╔═╡ aa0d9882-4093-11eb-1038-f3aa407bb8e9
md"# Airfoil Shape Design"

# ╔═╡ a4f74210-4093-11eb-10b3-c3ec247080d3
md"## Surface Optimization"

# ╔═╡ b618f06e-4093-11eb-17af-f918c12aa8da
function optimize_CST(αs, α, n_upper :: Integer, le = 0.)
    airfoil = (Foil ∘ kulfan_CST)(αs[1:n_upper], αs[n_upper+1:end], (0., 0.), le, 80)
    uniform = Uniform2D(1.0, α)
    cl = solve_case(airfoil, uniform, 60)
end

# ╔═╡ bc6f40f0-4093-11eb-30a1-17ffd1f0d74b
begin
	dzs = (0., 0.)
	α_u0 = [0.1, 0.1, 0.1, 0.1, 0.1]
	α_l0 = [-0.2, -0.1, -0.1, -0.1]
	α_ul0 = [α_u0; α_l0]

	n_upper = length(α_u0)  # Number of upper surface variables
	le = 0.0                # Leading edge modification
	α = 0.                  # Angle of attack
end

# ╔═╡ c5f9ba60-4093-11eb-0ea9-d5a3bf495234
CST_cl = optimize_CST(α_ul0, α, n_upper, le)

# ╔═╡ d7d91370-4093-11eb-0be4-01686a5cec7d
begin
	CST_test = kulfan_CST(α_u0, α_l0, dzs, le, 80)
	plot(CST_test[:,1], CST_test[:,2], 
		 label = "Initial Surface", aspectratio = 1)
end

# ╔═╡ 57bc7c80-4094-11eb-3e8b-1b496ce43efb
cl0 = 1.2               # Target lift coefficient

# ╔═╡ e8a268a0-4093-11eb-02d8-07eee9ba7d15
l_bound = [ fill(1e-12, length(α_u0)); fill(-Inf, length(α_l0)) ]

# ╔═╡ 2b8228e0-4094-11eb-174e-7b1f2d3fa6ad
u_bound = [ fill(Inf, length(α_u0)); fill(5e-1, length(α_l0)) ]

# ╔═╡ ec6e16a2-4093-11eb-1904-23cad91aff99
# resi_CST = optimize(x -> abs(optimize_CST(x, α, n_upper, le) - cl0), 
					# l_bound, u_bound,			# Bounds
					# α_ul0,						# Initial value
					# Fminbox(GradientDescent()),
					# autodiff = :forward,
					# Optim.Options(show_trace = true)
					# )

# ╔═╡ 0afe7d80-4094-11eb-0830-1d3a7e58c1db
begin
		CST_opt = kulfan_CST(resi_CST.minimizer[1:n_upper], resi_CST.minimizer[n_upper+1:end], dzs, le, 80)

	fig = plot(aspectratio = 1, dpi = 300)
	plot!(CST_test[:,1], CST_test[:,2], 
			label = "Initial Surface")
	plot!(CST_opt[:,1], CST_opt[:,2], 
			label = "Optimal Surface")
end

# ╔═╡ 8ed3fdc0-4093-11eb-1743-332b9934cbf8
## Plot


## Optimization



## Plot

# savefig(fig, "_research/tmp/CST_surface.pdf")
                
## Optimal test
# CST_cl = optimize_CST(resi_CST.minimizer, α, n_upper)
# println("CST Cl: $CST_cl")

# ## Camber-thickness optimization
# #============================================#

# # Camber optimization
# function optimize_camber_CST(αs, α, num_cam)
#     airfoil = (Foil ∘ camber_CST)(αs[1:num_cam], αs[num_cam+1:end], (0., 0.), 0., 80)
#     uniform = Uniform2D(1.0, α)
#     cl = solve_case(airfoil, uniform, 80)
# end

# ## Test
# α_c0 = zeros(6)
# α_t0 = fill(0.4, 6)
# α_ct0 = [α_c0; α_t0]

# num_cam = length(α_c0)  # Number of camber variables
# α = 0.                  # Angle of attack
# cl0 = 1.2               # Target lift coefficient
# camber_cl = optimize_camber_CST(α_ct0, α, num_cam)
# println("Camber Cl: $camber_cl")

# ## Plot
# camber_test = camber_CST(α_ct0[1:num_cam], α_ct0[num_cam+1:end], (0., 0.), 0., 80)
# camthick_test = foil_camthick(camber_test)
# plot(camber_test[:,1], camber_test[:,2], 
# 	label = "Initial Surface", aspectratio = 1)	
# # plot(camthick_test[:,1], camthick_test[:,2], 
# # 	label = "Initial Camber", aspectratio = 1)
# # plot!(camthick_test[:,1], camthick_test[:,3], 
# # 	label = "Initial Thickness", aspectratio = 1)

# ## Optimization
# l_bound = [ fill(-1e-12, length(α_c0)); [0.1, 0.12, 0.12, 0.05, 0.05, 0.05] ]
# u_bound = fill(Inf, length(α_ct0))

# resi_cam = optimize(x -> abs(optimize_camber_CST(x, α, num_cam) - cl0), 
# 					l_bound, u_bound,				# Bounds
# 					α_ct0,								# Initial value
# 					Fminbox(GradientDescent()),
# 					autodiff = :forward,
#                 	Optim.Options(
# 								#   extended_trace = true,
# 								  show_trace = true
# 								 )
# 					)

# ## Plot
# camber_opt = camber_CST(resi_cam.minimizer[1:num_cam], resi_cam.minimizer[num_cam+1:end], (0., 0.), 0., 80)

# surf_fig = plot(aspectratio = 1, dpi = 300)
# plot!(camber_test[:,1], camber_test[:,2], 
# 	label = "Initial Surface", aspectratio = 1)
# plot!(camber_opt[:,1], camber_opt[:,2], 
# 	label = "Optimal Surface for α = $α, Cl = $cl0")
# savefig(surf_fig, "_research/tmp/surface.pdf")

# ##
# camthick_opt = foil_camthick(camber_opt)
# camthick_fig = plot(aspectratio = 1, dpi = 300)
# plot!(camthick_test[:,1], camthick_test[:,2], 
# 	label = "Initial Camber", aspectratio = 1)
# plot!(camthick_test[:,1], camthick_test[:,3], 
# 	label = "Initial Thickness", aspectratio = 1)
# plot!(camthick_opt[:,1], camthick_opt[:,2], 
# 	label = "Optimal Camber for α = $α, Cl = $cl0", aspectratio = 1)
# plot!(camthick_opt[:,1], camthick_opt[:,3], 
# 	label = "Optimal Thickness for α = $α, Cl = $cl0", aspectratio = 1)
# savefig(camthick_fig, "_research/tmp/camthick.pdf")
                
# ## Optimal test
# camber_cl = optimize_camber_CST(resi_cam.minimizer, α, num_cam)
# println("Camber Cl: $camber_cl")

# ╔═╡ Cell order:
# ╠═aa0d9882-4093-11eb-1038-f3aa407bb8e9
# ╠═9c8e3e2e-4093-11eb-3e43-85d1fead00e5
# ╟─a4f74210-4093-11eb-10b3-c3ec247080d3
# ╠═b618f06e-4093-11eb-17af-f918c12aa8da
# ╠═bc6f40f0-4093-11eb-30a1-17ffd1f0d74b
# ╠═c5f9ba60-4093-11eb-0ea9-d5a3bf495234
# ╠═d7d91370-4093-11eb-0be4-01686a5cec7d
# ╠═57bc7c80-4094-11eb-3e8b-1b496ce43efb
# ╠═e8a268a0-4093-11eb-02d8-07eee9ba7d15
# ╠═2b8228e0-4094-11eb-174e-7b1f2d3fa6ad
# ╠═ec6e16a2-4093-11eb-1904-23cad91aff99
# ╠═0afe7d80-4094-11eb-0830-1d3a7e58c1db
# ╟─8ed3fdc0-4093-11eb-1743-332b9934cbf8
