### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d72bd840-3305-11eb-362b-0f6a7a774a13
begin
	using PlutoUI
	using StaticArrays
	using Rotations
	using Revise
	using AeroMDAO
end

# ╔═╡ bc905010-32f6-11eb-0557-91743b6fe9e3
using Plots; gr();

# ╔═╡ 23137140-3091-11eb-3f8a-2b8bd4b426b5
md"""
## Fixed Wing-VTOL UAV Sizing
"""

# ╔═╡ df735200-3263-11eb-3c15-c9f42ce9fbef
md"""
Tail sizing functions:
"""

# ╔═╡ 5fb5d550-3093-11eb-0691-09050e6939eb
md"""
```math
\begin{equation}
S_{HT} = \frac{C_{HT} \bar{c}_w S_w}{l_{HT}}
\end{equation}
```
"""

# ╔═╡ 304ebdb0-3091-11eb-0603-dd2131ef1226
horizontal_tail_area(V_h, mac_w, S_w, l_h) = V_h * mac_w * S_w / l_h;

# ╔═╡ 4da4b550-32ed-11eb-1709-1903b01b356c
md"""
```math
\begin{equation}
S_{VT} = \frac{C_{VT} b_w S_w}{2l_{VT}}
\end{equation}
```
"""

# ╔═╡ 4ba38e1e-32ed-11eb-2ed4-154da8d4afe8
vertical_tail_area(V_v, b_w, S_w, l_v) = V_v * b_w * S_w / 2l_v;

# ╔═╡ 6d9bd820-32ed-11eb-29c0-2b1f35984e1b
md"""
```math
\begin{equation}
b_{VT} = \frac{2 S_{VT}}{c_{HT}\left(1 + \frac{1}{\lambda_{VT}}\right)}
\end{equation}
```
"""

# ╔═╡ 4ba47880-32ed-11eb-206a-033b74ac2bc5
vertical_tail_span(S_v, c_h, λ_v) = 2 * S_v / c_h / (1 + 1. / λ_v);

# ╔═╡ 4ba53bd0-32ed-11eb-25b2-89385a3cecdc
htail_arm(x_VT, b_v, Λ_LE_v, c_h, x_CG) = x_VT + b_v * tan(Λ_LE_v) + 0.25 * c_h - x_CG;

# ╔═╡ 4bb0d490-32ed-11eb-1fe6-e72b49f18433
vtail_arm(x_mac_v, mac_v) = x_mac_v + 0.25 * mac_v;

# ╔═╡ f5ae15b2-3091-11eb-3c08-53653717475a
begin
	rect_chord(area, span) = area / span
	mean_aerodynamic_chord(c_r, λ) = 2/3 * c_r * (1 + λ + λ^2) / (1 + λ)
	aspect_ratio(span, area) = span^2 / area
end;

# ╔═╡ d5b5ec00-3092-11eb-2d10-87fcea8784a9
x_mac(x_cg, S_v, c_r_v, c_t_v) = x_cg + S_v * (c_r_v + 2 * c_t_v) / 3(c_r_v + c_t_v);

# ╔═╡ c2940800-32db-11eb-3dd7-6de7abd0d51c
function tail_arms(S_h, S_v, c_h, λ_v, Λ_LE_v, x_VT_mac, x_CG)
	b_v = vertical_tail_span(S_v, c_h, λ_v)				# Span
	c_r_v = c_h / λ_v 									# Root chord
	mac_v = mean_aerodynamic_chord(c_r_v, λ_v)			# MAC
	l_h = htail_arm(x_VT_mac, b_v, Λ_LE_v, c_h, x_CG)	# HTail arm
	l_v = vtail_arm(x_VT_mac, mac_v)					# VTail arm
	
	l_h, l_v
end;

# ╔═╡ 8cfe52f0-3267-11eb-2b11-23458f7f2966
md"""
Wing Parameters
"""

# ╔═╡ 6ca96a60-3304-11eb-150f-71f4ce3a2918
begin
	wing_foil = naca4((4,4,1,2))
	wing_num_secs = 3
	wing_foils = Foil.(Point2D.(first.(wing_foil), last.(wing_foil)) for i ∈ 1:wing_num_secs)
	
	wing_right = HalfWing(wing_foils,		# Foils
						 [0.2, 0.2, 0.1], 	# Chords
						 [0., 0., 0.],		# Twists
						 [1.705 / 2, 0.1], # Spans
						 [0., 60.], 		# Dihedrals
						 [0., 60.]) 		# Sweeps
	wing = Wing(wing_right, wing_right)
	info(wing)
end

# ╔═╡ 7f48f152-3263-11eb-36d1-cfd2a24b2cd9
begin
	x_w = 0.35
	Λ_LE_w = wing.right.sweeps[1]
	Λ_TE_w = deg2rad(0)
	c_r_w = wing.right.chords[1]
	c_t_w = wing.right.chords[2]
	λ_w =  taper_ratio(c_r_w, c_t_w)
	mac_w = mean_aerodynamic_chord(c_r_w, λ_w)
	b_w = span(wing)
	S_w = projected_area(wing)
end;

# ╔═╡ f11e9b60-3266-11eb-0c1f-fd0eaea6cc5d
V_H = 0.55

# ╔═╡ 8d97a350-30d3-11eb-307f-cd48a69b8300
V_V = 0.028

# ╔═╡ 98a56ace-3267-11eb-3627-ed82836734c2
md"""
Propeller parameters
"""

# ╔═╡ 3d4d2870-3263-11eb-3e80-1f3f5db82243
begin
	D_VTOL_prop = 0.36
	D_FW_prop = 0.16
	c_prop = D_VTOL_prop * 0.15
end;

# ╔═╡ 0fdf3ec2-32cf-11eb-399f-cb4281c72cc8
# Horizontal tail span
b_h = D_VTOL_prop + D_FW_prop

# ╔═╡ 322ba6c0-32d0-11eb-1557-dd4a0c9510b1
begin
	# Propeller locations
	x_PF = x_w + b_h / 2 * tan(Λ_LE_w) - (D_VTOL_prop / 2 + c_prop) / cos(Λ_LE_w)
	x_PR = x_w + c_r_w + b_h / 2 * tan(Λ_TE_w) + (D_VTOL_prop / 2 + c_prop) / cos(Λ_TE_w)
end;

# ╔═╡ 4d9a5b40-32d0-11eb-3977-c3a6511f65d1
begin
	# CG and vertical tail locations
	x_CG = 1/2 * (x_PR - x_PF)
	x_VT = x_PR + D_VTOL_prop / 2 + c_prop
end;

# ╔═╡ 40ba59f2-32eb-11eb-2b41-91c19ca0cf8c
function tail_sizing(num_iter, l_h, l_v, V_V, V_H, S_w, b_w, mac_w, λ_v, Λ_LE_v)
	# Initialisations
	S_h0, S_v0 = 1, 1
	S_hs = [ S_h0; zeros(num_iter - 1) ]
	S_vs = [ S_v0; zeros(num_iter - 1) ]
	error_S_h = [ 100; zeros(num_iter-1) ]
	error_S_v = [ 100; zeros(num_iter-1) ]
	
	for i ∈ 2:num_iter
		S_h = horizontal_tail_area(V_H, mac_w, S_w, l_h)
		S_v = vertical_tail_area(V_V, b_w, S_w, l_v)

		c_h = rect_chord(S_h, b_h)

		l_h, l_v = tail_arms(S_h, S_v, c_h, λ_v, Λ_LE_v, x_VT, x_CG)

		diffs = abs.([ (S_h - S_hs[i-1])/S_hs[i-1], (S_v - S_vs[i-1])/S_vs[i-1] ] .* 100)
		
		error_S_h[i] = diffs[1]
		error_S_v[i] = diffs[2]
		S_hs[i] = S_h
		S_vs[i] = S_v
	end
	
	S_hs, S_vs, l_h, l_v, error_S_h, error_S_v
end;

# ╔═╡ 34a80ae0-32f0-11eb-1dc3-07cae4d044ce
md"""
Iterations:
"""

# ╔═╡ fdc41530-32ce-11eb-22c0-a9bd5a40611d
begin
	# Vertical tail parameters
	λ_v = 0.5
	Λ_LE_v = 30.
end;

# ╔═╡ ab52a460-349a-11eb-1028-27f26f3c9606
max_iter = 15;

# ╔═╡ 82248a1e-33fc-11eb-3cc1-9b4ae3379bef
begin
	num_slider = @bind num NumberField(1:max_iter, default = 6)
	md"""
	Number of iterations: $(num_slider)
	"""
end

# ╔═╡ f926e9c0-32f7-11eb-1e53-3ba8e9f52056
S_hs, S_vs, l_h, l_v, error_h, error_v = tail_sizing(num, x_VT - x_CG, x_VT - x_CG, V_V, V_H, S_w, b_w, mac_w, λ_v, deg2rad(Λ_LE_v));

# ╔═╡ 8503ef50-332b-11eb-2953-03dc7826b417
begin
	c_h = rect_chord(S_hs[end], b_h)
	b_v = vertical_tail_span(S_vs[end], c_h, λ_v)
	c_r_v = c_h / λ_v
	mac_v = mean_aerodynamic_chord(c_r_v, λ_v)
end;

# ╔═╡ c1522370-332e-11eb-258f-a96268d86a87
begin
	htail_foil = naca4((0,0,1,2))
	htail_foils = Foil.(Point2D.(first.(htail_foil), last.(htail_foil)) for i ∈ 1:2)
	
	htail_right = HalfWing(htail_foils,
						  [c_h, c_h], 	# Chords
						  [0., 0.],		# Twists
						  [b_h / 2],  	# Spans
						  [0.],			# Dihedrals
						  [0.])			# Sweeps

	htail = Wing(htail_right, htail_right)
	info(htail)
end

# ╔═╡ 7dcd3140-3328-11eb-2e98-fbc13dbfcfa4
begin
	vtail_foil = naca4((0,0,0,9))
	vtail_foils = Foil.(Point2D.(first.(vtail_foil), last.(vtail_foil))  for i ∈ 1:2)
	
	vtail_left = HalfWing(vtail_foils,
						 [c_r_v, λ_v * c_r_v], 	# Chords
						 [0., 0.],					# Twists				 
						 [b_v],  					# Spans
						 [0.],						# Dihedrals
						 [Λ_LE_v])					# Sweeps
						 
	vtail_right = vtail_left
	info(vtail_left)
end

# ╔═╡ c06dff70-32f6-11eb-2a38-8b82b14cebee
begin
	plot(1:length(error_h), error_h, label = "Horizontal Tail Error")
	plot!(1:length(error_v), error_v, label = "Vertical Tail Error")
end

# ╔═╡ c4943450-34c9-11eb-07f6-ab75cfb6de6e
begin
	φ_s = @bind φ Slider(0:1e-2:90, default = 15)
	ψ_s = @bind ψ Slider(0:1e-2:90, default = 30)
	z_s = @bind z_limit Slider(0:1e-2:3, default = 1.)
	md"""
	Horizontal: $(φ_s)
	Vertical: $(ψ_s)
	
	*z*-scale: $(z_s)
	"""
end

# ╔═╡ 885f6532-33c1-11eb-1b81-e1a81e959643
begin
	leading, trailing = tupvector.(wing_bounds(wing))
	leading_h, trailing_h = tupvector.(wing_bounds(htail))
	leading_vr, trailing_vr = tupvector.(wing_bounds(vtail_right))
	leading_vl, trailing_vl = tupvector.(wing_bounds(vtail_left))
end;

# ╔═╡ 8f44ca92-33a6-11eb-315f-5b87f51ef53c
begin
	wing_coords = [ (x_w, 0, 0) .+ coords for coords in [ leading; trailing[end:-1:1]; leading[1] ] ][:]
	htail_coords = [ (x_CG + l_h - 0.25 * c_h, 0, b_v) .+ coords for coords in [ leading_h; trailing_h[end:-1:1]; leading_h[1] ] ][:]
	vtail1_coords = tupvector(RotX(π/2) * SVector(coords...) .+ SVector(x_VT, -b_h / 2, 0) for coords in [ leading_vl; trailing_vl[end:-1:1]; leading_vl[1] ])[:]
	vtail2_coords = tupvector(RotX(π/2) * SVector(coords...) .+ SVector(x_VT, b_h / 2, 0) for coords in [ leading_vr; trailing_vr[end:-1:1]; leading_vr[1] ])[:]
end;

# ╔═╡ 2c39cf70-33b7-11eb-0de3-236d16a6f5ea
circle3D(r) = [ (r*cos(θ), r*sin(θ), 0) for θ in 0:1e-2:2π ];

# ╔═╡ 4bc3ffa0-33b7-11eb-0525-659b21f803eb
begin
	circ3D = circle3D(D_VTOL_prop / 2)
	circ3D_fw = circle3D(D_FW_prop / 2)
	prop3D_rear_right = [ (x_PR, (b_h)/ 2, 0) .+ coords for coords in circ3D ]
	prop3D_rear_left = [ (x_PR, -(b_h) / 2, 0) .+ coords for coords in circ3D ]
	prop3D_front_left = [ (x_PF, -(b_h) / 2, 0) .+ coords for coords in circ3D ]
	prop3D_front_right = [ (x_PF, (b_h) / 2, 0) .+ coords for coords in circ3D ]
	prop3D_fw = tupvector((x_w + c_r_w * 1.1, 0, 0) .+ RotY(π/2) * SVector(coords...) for coords in circ3D_fw)
end;

# ╔═╡ 792a3af0-3479-11eb-37af-dda9a6268480
begin
	cuck = plot(camera = (φ, ψ), zlim=(-z_limit, z_limit), aspect_ratio = 1)
	plot!(wing_coords, label = "Wing")
	plot!(htail_coords, label = "Horizontal Tail")
	plot!(vtail1_coords, label = "Vertical Tail 1")
	plot!(vtail2_coords, label = "Vertical Tail 2")
	plot!(prop3D_rear_right, label = "Prop Rear Right")
	plot!(prop3D_rear_left, label = "Prop Rear Left")
	plot!(prop3D_front_left, label = "Prop Front Left")
	plot!(prop3D_front_right, label = "Prop Front Right")
	plot!(prop3D_fw, label = "Prop Fixed-Wing")
end

# ╔═╡ d31ac9a0-34c2-11eb-1dd0-0be54d9e8e1d
top_view = plot(cuck, camera = (0, 90))

# ╔═╡ 8523e000-34c3-11eb-0591-d5669fed462b
side_view = plot(cuck, camera = (0, 0))

# ╔═╡ Cell order:
# ╟─23137140-3091-11eb-3f8a-2b8bd4b426b5
# ╠═d72bd840-3305-11eb-362b-0f6a7a774a13
# ╟─df735200-3263-11eb-3c15-c9f42ce9fbef
# ╟─5fb5d550-3093-11eb-0691-09050e6939eb
# ╟─304ebdb0-3091-11eb-0603-dd2131ef1226
# ╟─4da4b550-32ed-11eb-1709-1903b01b356c
# ╟─4ba38e1e-32ed-11eb-2ed4-154da8d4afe8
# ╟─6d9bd820-32ed-11eb-29c0-2b1f35984e1b
# ╟─4ba47880-32ed-11eb-206a-033b74ac2bc5
# ╠═4ba53bd0-32ed-11eb-25b2-89385a3cecdc
# ╠═4bb0d490-32ed-11eb-1fe6-e72b49f18433
# ╠═f5ae15b2-3091-11eb-3c08-53653717475a
# ╠═d5b5ec00-3092-11eb-2d10-87fcea8784a9
# ╠═c2940800-32db-11eb-3dd7-6de7abd0d51c
# ╠═40ba59f2-32eb-11eb-2b41-91c19ca0cf8c
# ╟─8cfe52f0-3267-11eb-2b11-23458f7f2966
# ╠═6ca96a60-3304-11eb-150f-71f4ce3a2918
# ╠═7f48f152-3263-11eb-36d1-cfd2a24b2cd9
# ╠═f11e9b60-3266-11eb-0c1f-fd0eaea6cc5d
# ╠═8d97a350-30d3-11eb-307f-cd48a69b8300
# ╟─98a56ace-3267-11eb-3627-ed82836734c2
# ╠═3d4d2870-3263-11eb-3e80-1f3f5db82243
# ╟─0fdf3ec2-32cf-11eb-399f-cb4281c72cc8
# ╠═322ba6c0-32d0-11eb-1557-dd4a0c9510b1
# ╠═4d9a5b40-32d0-11eb-3977-c3a6511f65d1
# ╟─34a80ae0-32f0-11eb-1dc3-07cae4d044ce
# ╠═fdc41530-32ce-11eb-22c0-a9bd5a40611d
# ╠═f926e9c0-32f7-11eb-1e53-3ba8e9f52056
# ╠═8503ef50-332b-11eb-2953-03dc7826b417
# ╟─c1522370-332e-11eb-258f-a96268d86a87
# ╟─7dcd3140-3328-11eb-2e98-fbc13dbfcfa4
# ╠═bc905010-32f6-11eb-0557-91743b6fe9e3
# ╠═ab52a460-349a-11eb-1028-27f26f3c9606
# ╟─82248a1e-33fc-11eb-3cc1-9b4ae3379bef
# ╟─c06dff70-32f6-11eb-2a38-8b82b14cebee
# ╟─c4943450-34c9-11eb-07f6-ab75cfb6de6e
# ╟─792a3af0-3479-11eb-37af-dda9a6268480
# ╟─d31ac9a0-34c2-11eb-1dd0-0be54d9e8e1d
# ╟─8523e000-34c3-11eb-0591-d5669fed462b
# ╟─8f44ca92-33a6-11eb-315f-5b87f51ef53c
# ╟─4bc3ffa0-33b7-11eb-0525-659b21f803eb
# ╟─885f6532-33c1-11eb-1b81-e1a81e959643
# ╟─2c39cf70-33b7-11eb-0de3-236d16a6f5ea
