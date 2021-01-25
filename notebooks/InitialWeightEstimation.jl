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

# ╔═╡ 6c04f612-4083-11eb-017a-0f6e29496ddd
using Plots

# ╔═╡ caca53a0-40f3-11eb-211d-45d2e919cda0
using PlutoUI

# ╔═╡ 5c4d5f00-4083-11eb-31ad-d70d906e4825
md"""
# Initial Weight Estimation
![](https://godot-bloggy.xyz/WeightEstimation.svg)
"""

# ╔═╡ 2ede8630-4087-11eb-0f6b-4ddfbf717726
md"""## Maximum Takeoff Weight
![](https://godot-bloggy.xyz/XDSMTakeoff.svg)

```math
W_0 = \frac{W_\text{payload} + W_\text{crew}}{1 - W_{f_N} / W_0 - W_e / W_0}
```
"""

# ╔═╡ 0bd8d880-4103-11eb-01aa-6f908fa8cd9e
hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));

# ╔═╡ c2c891f0-4083-11eb-1164-f78e70d8cd88
"""
	maximum_takeoff_weight(payload_weight, crew_weight, fuel_weight_frac, empty_weight_frac)

Computes the maximum takeoff weight.
"""
maximum_takeoff_weight(WPL, Wcrew, WfWTO, WeWTO) = (WPL + Wcrew)/(1 - WfWTO - WeWTO)

# ╔═╡ 49a58902-4087-11eb-27ae-61c820955b7a
md"Gravitational acceleration constant"

# ╔═╡ 82e6a100-4086-11eb-0ee0-197bf26577e1
g = 9.81;

# ╔═╡ 68635df2-4086-11eb-2c0c-fd451c7480b9
md"### Takeoff"

# ╔═╡ bc194170-4082-11eb-3f1b-27d8e7490b12
takeoffWF = 0.97;

# ╔═╡ 706b07f0-4086-11eb-358d-9f296a4bf423
md"### Climb"

# ╔═╡ 4df075f0-4083-11eb-0f64-c5fd6c41772f
climbWF = 0.985;

# ╔═╡ 4f886920-4085-11eb-3720-45a3b9fa0a37
md"### Cruise

$$WF_{\mathrm{cruise}} = \exp\left(-\frac{R \times SFC}{V \times L/D}\right)$$"

# ╔═╡ 9682d72e-4084-11eb-269d-0bc3efc1c8e2
cruise_weight_fraction(range, SFC, V, L_D) = exp(-range * SFC / (V * L_D) )

# ╔═╡ 8f3f4b10-4085-11eb-36f7-0de9a630e3e3
begin
	# At 35000 ft speed of sound = 295 m/s
	M = 0.84 					# Cruise speed in Mach
	V = M * 295       			# Cruise speed in m/s
	cruise_SFC = 0.5/3600 		# TSFC at cruise in 1/secs
	R1 = 2800 * 1000    		# Range of cruise segment 1
	LD_max = 16.                # Maximum L/D ratio
	LD_cruise = LD_max * 0.866  # L/D ratio at cruise
end;

# ╔═╡ 36093a82-4083-11eb-346a-23528f9cd4d5
cruise1WF = cruise_weight_fraction(R1, cruise_SFC, V, LD_cruise)

# ╔═╡ 54b79ec0-4085-11eb-3386-41981e608e45
md"### Loiter

$$WF_{\mathrm{loiter}} = \exp\left(-\frac{E \times SFC}{L/D}\right)$$"

# ╔═╡ e956f3a0-5570-11eb-2311-59cedf436b9a
md"""Call your function `loiter_weight_fraction`!"""

# ╔═╡ bcb4fa50-4084-11eb-3474-618341f94e86
loiter_weight_fraction(endurance, SFC, L_D) = exp(-endurance * SFC / L_D)

# ╔═╡ 3afb38b0-4086-11eb-31c0-cd32f6feeed6
begin
	E1 = 3 * 3600
	loiter_SFC = 0.4/3600  
end;

# ╔═╡ 3b336710-4083-11eb-3990-8ff6fcd8c645
loiter1WF = 1.

# ╔═╡ 2db14640-4086-11eb-34ab-7982af7dd0b5
md"Range of cruise segment 2:"

# ╔═╡ 27ee5c20-4086-11eb-25e6-3ff461fb24dd
R2 = 2800 * 1000;

# ╔═╡ f8b34110-4084-11eb-2b95-2fb540813285
cruise2WF = cruise_weight_fraction(R2, cruise_SFC, V, LD_cruise)

# ╔═╡ 191983f0-4086-11eb-0575-7d769008814e
md"Endurance of loiter segment 2 in seconds:"

# ╔═╡ 1197d2d0-4086-11eb-2134-4bf31d078ddc
E2 = 0.33 * 3600

# ╔═╡ 07ecf590-4085-11eb-1181-633ab1ab1e4e
loiter2WF = loiter_weight_fraction(E2, loiter_SFC, LD_max)

# ╔═╡ 0ec73130-4087-11eb-064c-f5695f578c69
md"### Landing"

# ╔═╡ 4950af0e-4083-11eb-2ec1-1dd821c489b4
landingWF = 0.995;

# ╔═╡ 158ee120-4087-11eb-188b-f13671908d3f
md"""## Fuel Weight Fractions
![](https://godot-bloggy.xyz/XDSMFuel.svg)
```math
W_{f_N} / W_0 = a\left(1 - \prod_{i = 1}^{N}\frac{W_{f_i}}{W_{f_{i-1}}}\right)
```
"""

# ╔═╡ edc9887e-4102-11eb-1452-7336ed39e616
hint(md"""
	You can evaluate the product of an array using `prod()`!
	""")

# ╔═╡ 08eb3da0-4083-11eb-2eac-45e93f3e99f4
fuel_weight_fraction(fracs, a = 1.00) = a * (1 - prod(fracs))

# ╔═╡ 2816e82e-4085-11eb-33d7-93dac498b4ff
md"Create list of fuel fractions"

# ╔═╡ 0e95bdc0-4083-11eb-267c-07ddf6161aa6
FFs = [takeoffWF, climbWF, cruise1WF, loiter1WF, cruise2WF, loiter2WF, landingWF]

# ╔═╡ 13bc67e2-4083-11eb-1b48-d5fa25a4cd95
WfWTO = fuel_weight_fraction(FFs)

# ╔═╡ 3bd51360-4085-11eb-240e-e189c1d79c9c
md"""## Empty Weight Fraction
![](https://godot-bloggy.xyz/XDSMEmpty.svg)
"""

# ╔═╡ ba7a8fb0-40fd-11eb-1d5e-7bd861fae148
md"""
Raymer: $W_e / W_0 = A W_0^B$
"""

# ╔═╡ c6bbec32-4083-11eb-32a6-93488d9abc07
function empty_weight_raymer(WTO, A, B)
	WeWTO = A * WTO^B   # Raymer's regression
	We = WeWTO * WTO     # Empty weight calculation
	return WeWTO
end

# ╔═╡ df1745c0-40fd-11eb-21ed-3b2ebf4f1dbd
md"""
Roskam: $W_e / W_0 = 10^{(\log W_0 - A) / B}$
"""

# ╔═╡ cbe0c190-4083-11eb-020e-bbf9f4f7d67a
function empty_weight_roskam(WTO, A, B)
	logWe = (log10(WTO) - A) / B    # Roskam's regression
	We = 10^logWe
	WeWTO = We/WTO                     # Empty weight fraction calculation
	return We, WeWTO
end

# ╔═╡ bb6c4652-4087-11eb-0818-2515b67bccc0
md"Regression coefficients"

# ╔═╡ d9628ce0-4083-11eb-2d1d-c951a4c14425
begin
	A = 0.88
	B = -0.07
end;

# ╔═╡ 33ab3060-40fe-11eb-00d9-071f7ca836e9
md"""
## Iteration
![](https://godot-bloggy.xyz/XDSMTakeoffEmpty.svg)
"""

# ╔═╡ 089674a0-4088-11eb-1534-ef06a695815f
function compute_mtow(W_PL, W_crew, A, B, num_iters = 20, err = 1.0, tol = 1e-12)
	WTO = W_PL + W_crew 	# Initial value (guess)
	WTOs = [WTO]
	errors = []
	for i in 1:num_iters
		WeWTO = empty_weight_raymer(WTO, A, B)
		newWTO = maximum_takeoff_weight(W_PL, W_crew, WfWTO, WeWTO)
		error = abs(newWTO-WTO)/WTO
		WTO = newWTO
		push!(WTOs, WTO)
		push!(errors, error)
		error < tol ? break : continue
	end
	
	return WTOs, errors
end

# ╔═╡ 89e9ac40-4086-11eb-3ca3-f1352e475899
begin
	WPL = 4500g                           # Payload weight in N
	Wcrew = 360g                          # Crew weight in N
end;

# ╔═╡ bf757a70-40f3-11eb-3948-31c68d44746d
max_iter = 20

# ╔═╡ b2820c70-40f3-11eb-0600-037e8d844797
begin
	num_slider = @bind num NumberField(1:max_iter, default = 10)
	md"""
	Number of iterations: $(num_slider)
	"""
end

# ╔═╡ 05618530-4084-11eb-325c-1fe833f089dc
WTOs, errors = compute_mtow(WPL, Wcrew, A, B, num);

# ╔═╡ 58da9050-408c-11eb-310e-6106a7b11edd
begin
	plot1 = plot(WTOs, label = :none, ylabel = "MTOW", xlabel = "Iterations")
	plot2 = plot(errors, label = :none, ylabel = "Error", xlabel = "Iterations")
	plot(plot1, plot2, layout = (2,1))
end

# ╔═╡ Cell order:
# ╟─5c4d5f00-4083-11eb-31ad-d70d906e4825
# ╟─2ede8630-4087-11eb-0f6b-4ddfbf717726
# ╟─0bd8d880-4103-11eb-01aa-6f908fa8cd9e
# ╠═c2c891f0-4083-11eb-1164-f78e70d8cd88
# ╟─49a58902-4087-11eb-27ae-61c820955b7a
# ╠═82e6a100-4086-11eb-0ee0-197bf26577e1
# ╟─68635df2-4086-11eb-2c0c-fd451c7480b9
# ╠═bc194170-4082-11eb-3f1b-27d8e7490b12
# ╟─706b07f0-4086-11eb-358d-9f296a4bf423
# ╠═4df075f0-4083-11eb-0f64-c5fd6c41772f
# ╟─4f886920-4085-11eb-3720-45a3b9fa0a37
# ╟─9682d72e-4084-11eb-269d-0bc3efc1c8e2
# ╠═8f3f4b10-4085-11eb-36f7-0de9a630e3e3
# ╟─36093a82-4083-11eb-346a-23528f9cd4d5
# ╟─54b79ec0-4085-11eb-3386-41981e608e45
# ╟─e956f3a0-5570-11eb-2311-59cedf436b9a
# ╟─bcb4fa50-4084-11eb-3474-618341f94e86
# ╟─3afb38b0-4086-11eb-31c0-cd32f6feeed6
# ╠═3b336710-4083-11eb-3990-8ff6fcd8c645
# ╟─2db14640-4086-11eb-34ab-7982af7dd0b5
# ╠═27ee5c20-4086-11eb-25e6-3ff461fb24dd
# ╠═f8b34110-4084-11eb-2b95-2fb540813285
# ╟─191983f0-4086-11eb-0575-7d769008814e
# ╟─1197d2d0-4086-11eb-2134-4bf31d078ddc
# ╠═07ecf590-4085-11eb-1181-633ab1ab1e4e
# ╟─0ec73130-4087-11eb-064c-f5695f578c69
# ╠═4950af0e-4083-11eb-2ec1-1dd821c489b4
# ╟─158ee120-4087-11eb-188b-f13671908d3f
# ╟─edc9887e-4102-11eb-1452-7336ed39e616
# ╟─08eb3da0-4083-11eb-2eac-45e93f3e99f4
# ╟─2816e82e-4085-11eb-33d7-93dac498b4ff
# ╠═0e95bdc0-4083-11eb-267c-07ddf6161aa6
# ╠═13bc67e2-4083-11eb-1b48-d5fa25a4cd95
# ╟─3bd51360-4085-11eb-240e-e189c1d79c9c
# ╟─ba7a8fb0-40fd-11eb-1d5e-7bd861fae148
# ╟─c6bbec32-4083-11eb-32a6-93488d9abc07
# ╟─df1745c0-40fd-11eb-21ed-3b2ebf4f1dbd
# ╟─cbe0c190-4083-11eb-020e-bbf9f4f7d67a
# ╟─bb6c4652-4087-11eb-0818-2515b67bccc0
# ╠═d9628ce0-4083-11eb-2d1d-c951a4c14425
# ╟─33ab3060-40fe-11eb-00d9-071f7ca836e9
# ╟─089674a0-4088-11eb-1534-ef06a695815f
# ╠═89e9ac40-4086-11eb-3ca3-f1352e475899
# ╠═bf757a70-40f3-11eb-3948-31c68d44746d
# ╟─b2820c70-40f3-11eb-0600-037e8d844797
# ╠═05618530-4084-11eb-325c-1fe833f089dc
# ╟─6c04f612-4083-11eb-017a-0f6e29496ddd
# ╟─58da9050-408c-11eb-310e-6106a7b11edd
# ╟─caca53a0-40f3-11eb-211d-45d2e919cda0
