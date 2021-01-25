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

# ╔═╡ f7132010-404c-11eb-157e-0dfc3c36ead3
begin
	using PlutoUI
	using Revise
	using StaticArrays
	using AeroMDAO
	using Plots
	gr()
end;

# ╔═╡ e9306402-5ee0-11eb-1cbb-8dd6f2e08ad6
md"""
# Vortex Lattice Method
"""

# ╔═╡ 63bd2ac0-5f0c-11eb-0374-f7584e911a22
md"## Wing Definition"

# ╔═╡ 0a99e790-404d-11eb-2b93-f3d0acd94967
foil = naca4((2,4,1,2));

# ╔═╡ 0aa950e0-404d-11eb-3a8a-bfe40df19222
wing_right = HalfWing(Foil.(foil for i ∈ 1:3),	# Foils
                      [0.4, 0.2, 0.1], 			# Chords
					  [0., 2., 5.],				# Twists
					  [1.705 / 2, 0.1],			# Spans
					  [0., 60.],				# Dihedrals
					  [0., 30.]);				# Sweeps

# ╔═╡ 0aae0bd0-404d-11eb-1245-9186c5b9c21f
wing = Wing(wing_right, wing_right);

# ╔═╡ a1a13230-5ef0-11eb-3ee9-99f383b5cbd7
wing_span, wing_area, wing_mac, wing_AR = info(wing);

# ╔═╡ 0aaf4450-404d-11eb-31ba-757f8422c1ab
md"""
Parameter | Value
:-------- | -----:
Span ($m$)     | $wing_span
Planform Area ($m^2$) | $wing_area 
Mean Aerodynamic Chord ($m$) | $wing_mac
Aspect Ratio | $wing_AR
"""

# ╔═╡ 6bf6e140-5f0c-11eb-3cdc-4744f2861b67
md"## Analysis Setup"

# ╔═╡ 03f01a5e-5f0e-11eb-253a-135b7e33a9a4
md"For an analysis, you require ambient and conditions. In this case, you need the density $\rho$, a reference location $x_\text{ref}$ for calculating moments."

# ╔═╡ 0ac53d50-404d-11eb-0485-637991172933
begin
	ρ = 1.225
	ref = SVector(0.25 * mean_aerodynamic_chord(wing), 0., 0.)
end;

# ╔═╡ 436fe6c0-5f0e-11eb-38a7-696752b4eebf
md"""

A `Freestream` object requires the boundary conditions of the analysis, given by the following:

Parameter |  Value  | Description
:-------- | :-----: | ----------:
$V_\infty$|   $10$    | Freestream speed
$\alpha$  |    $5$    | Angle of attack (deg)
$\beta$   |    $5$    | Sideslip angle (deg)
$\Omega$  | $(0,0,0)$ | Rotation vector
"""

# ╔═╡ 3c590890-5f12-11eb-05f9-03636fc8085a
V, α, β = 10.0, 5.0, 0.0

# ╔═╡ c5631d00-5f0e-11eb-0878-f1f1cdf50b5f
begin
	Ω = SVector(0.0, 0.0, 0.0)
	uniform = Freestream(V, α, β, Ω)
end;

# ╔═╡ 0add0b10-404d-11eb-0b56-91f9e3431c70
nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 15, chord_num = 10);

# ╔═╡ 5889d8c0-5f0b-11eb-174f-d7a48f04186a
begin
	CL_nf, CDi_nf, CY_nf, Cl_nf, Cm_nf, Cn_nf, p_b_nf, q_b_nf, r_b_nf = nf_coeffs
	CL_ff, CDi_ff, CY_ff, Cl_ff, Cm_ff, Cn_ff, p_b_ff, q_b_ff, r_b_ff = ff_coeffs
end;

# ╔═╡ 81833050-404d-11eb-0dbd-7fc5e4d69d77
md"""
Coefficient | Nearfield Value | Farfield Value
:---------- | :-------------: |  -------------:
$C_L$     	|	$CL_nf        | $CL_ff
$C_{D_i}$ 	|	$CDi_nf       | $CDi_ff
$C_Y$     	|	$CY_nf 		  | $CY_ff
$C_l$       |	$Cl_nf 		  | $Cl_ff
$C_m$		|	$Cm_nf        | $Cm_ff
$C_n$		|	$Cn_nf        | $Cn_ff
$\bar{p}$	|	$p_b_nf       | $p_b_ff
$\bar{q}$	|	$q_b_nf       | $q_b_ff
$\bar{r}$	|	$r_b_nf       | $r_b_ff

"""

# ╔═╡ 783222d0-5f0c-11eb-3ba6-0533ee09b40c
md"## Visualization"

# ╔═╡ 6969d4f2-5ee4-11eb-1994-f9da9720641f
md"Seed 1"

# ╔═╡ 50efe150-5ee2-11eb-2c3a-494788215ba7
begin
	num_points = 50
	max_z = 0.1
	y = span(wing) / 2 - 0.01
	seed_1 = SVector.(fill(-0.1, num_points), fill(y, num_points), range(-max_z, stop = max_z, length = num_points))
end;

# ╔═╡ 63598e20-5ee4-11eb-20ed-81a6a711d0ec
md"Seed 2"

# ╔═╡ 598d6d40-5ee3-11eb-35e0-5725d603e357
begin
	span_points = 20
	init        = trailing_chopper(ifelse(β == 0, wing.right, wing), span_points) 
	dx, dy, dz  = 0, 0, 1e-3
	seed_2        = [ init .+ Ref([dx, dy, dz])  ; 
					  init .+ Ref([dx, dy, -dz]) ]
end;

# ╔═╡ 0b11fdc2-404d-11eb-269c-6904c57e599d
begin
	horseshoe_coords = plot_panels(horseshoe_panels[:])
	camber_coords = plot_panels(camber_panels[:])
	wing_coords = plot_surface(wing);
	streams = plot_streams(uniform, seed_2, horseshoes, Γs, 2, 100)
end;

# ╔═╡ 158416d0-40f7-11eb-34a9-c193e6456ff3
begin
	φ_s = @bind φ Slider(0:1e-2:90, default = 15)
	ψ_s = @bind ψ Slider(0:1e-2:90, default = 30)
	z_s = @bind z_limit Slider(0:1e-2:3, default = 1.)
	stream = @bind stream CheckBox()
	md"""
	Horizontal: $(φ_s)
	Vertical: $(ψ_s)
	
	*z*-scale: $(z_s)
	Streamlines: $stream
	"""
end

# ╔═╡ bcf98ed0-5f0c-11eb-1fac-1ff84a92a08b
begin
	plot(xaxis = "x", yaxis = "y", zaxis = "z",
			 aspect_ratio = 1, 
			 camera = (φ, ψ),
			 zlim = (-0.1, z_limit)
			)
		plot!.(wing_coords, color = :black, label = :none)
	plot!()
end

# ╔═╡ 0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
begin
	plot(xaxis = "x", yaxis = "y", zaxis = "z",
		 aspect_ratio = 1, 
		 camera = (φ, ψ),
		 zlim = (-0.1, z_limit)
		)
	plot!.(horseshoe_coords, color = :black, label = :none)
	if stream
		plot!.(streams, color = :green, label = :none)
	end
	plot!()
end

# ╔═╡ Cell order:
# ╟─e9306402-5ee0-11eb-1cbb-8dd6f2e08ad6
# ╟─63bd2ac0-5f0c-11eb-0374-f7584e911a22
# ╠═0a99e790-404d-11eb-2b93-f3d0acd94967
# ╠═0aa950e0-404d-11eb-3a8a-bfe40df19222
# ╠═0aae0bd0-404d-11eb-1245-9186c5b9c21f
# ╠═a1a13230-5ef0-11eb-3ee9-99f383b5cbd7
# ╠═0aaf4450-404d-11eb-31ba-757f8422c1ab
# ╟─bcf98ed0-5f0c-11eb-1fac-1ff84a92a08b
# ╟─6bf6e140-5f0c-11eb-3cdc-4744f2861b67
# ╟─03f01a5e-5f0e-11eb-253a-135b7e33a9a4
# ╠═0ac53d50-404d-11eb-0485-637991172933
# ╟─436fe6c0-5f0e-11eb-38a7-696752b4eebf
# ╠═3c590890-5f12-11eb-05f9-03636fc8085a
# ╠═c5631d00-5f0e-11eb-0878-f1f1cdf50b5f
# ╠═0add0b10-404d-11eb-0b56-91f9e3431c70
# ╠═5889d8c0-5f0b-11eb-174f-d7a48f04186a
# ╟─81833050-404d-11eb-0dbd-7fc5e4d69d77
# ╟─783222d0-5f0c-11eb-3ba6-0533ee09b40c
# ╟─6969d4f2-5ee4-11eb-1994-f9da9720641f
# ╠═50efe150-5ee2-11eb-2c3a-494788215ba7
# ╟─63598e20-5ee4-11eb-20ed-81a6a711d0ec
# ╠═598d6d40-5ee3-11eb-35e0-5725d603e357
# ╟─0b11fdc2-404d-11eb-269c-6904c57e599d
# ╟─158416d0-40f7-11eb-34a9-c193e6456ff3
# ╟─0b4346f0-404d-11eb-1af9-e9b0f2a7ee21
# ╟─f7132010-404c-11eb-157e-0dfc3c36ead3
