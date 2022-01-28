### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 75407f80-819a-11eb-3307-e5c7d2dcef46
using AeroMDAO

# ╔═╡ 7542a260-819a-11eb-36b0-4d6f6497e611
using Plots

# ╔═╡ 75665700-819a-11eb-04e3-23d133e8dc6c
using StaticArrays

# ╔═╡ 91040ac0-819a-11eb-2967-2909f8487f5b
using PlutoUI

# ╔═╡ 993a3ed2-819a-11eb-0464-ab84db58a8a2
md"""
# AeroMDAO -- Airfoil Analyses
**References**: 

1. Mark Drela. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
2. Joseph Katz and Allen Plotkin. _Low-Speed Aerodynamics, Second Edition_. Cambridge University Press, 2001.
"""

# ╔═╡ 526c4d80-f348-448a-a8f0-ff073bd96418
gr(
	grid = true, 	   # Disable the grid for a e s t h e t i c
	size = (800, 520), # Adjust the dimensions to suit your monitor
	dpi = 300     	   # Higher DPI = higher density
  )

# ╔═╡ 754bca22-819a-11eb-0039-59eb17946634
md"## Your First Airfoil"

# ╔═╡ 2060cc90-819f-11eb-1b1e-d1f61e2c878a
md"### NACA 4--digit Airfoils"

# ╔═╡ 754cdb90-819a-11eb-30f7-0d172a955548
md"Your can define a NACA-4 airfoil using the following function:"

# ╔═╡ 754dc5ee-819a-11eb-373d-f909ce6e2e3a
md"""
```julia
naca4(digits 		:: NTuple{4, Real}, # Digits, e.g. (2,4,1,2)
	  points = 40 	:: Integer; 		# Number of points
      sharp_trailing_edge = true)		# Sharp or blunt trailing edge
```
"""

# ╔═╡ 755455a2-819a-11eb-1384-a184f1a5ed93
digits = (4,4,1,2);

# ╔═╡ 7554f1e0-819a-11eb-06bb-6587881bb0a1
airfoil = naca4(digits, 60; sharp_trailing_edge = true)

# ╔═╡ 758a32b0-819a-11eb-385d-e75b23f41ac4
md"### Geometric Representations"

# ╔═╡ 758b6b30-819a-11eb-35c6-99ad31016f75
md"You can convert these coordinates into the camber-thickness representation by calling:

```julia
coordinates_to_camber_thickness(coords, 						# Coordinates
			  num_points = 60 :: Integer)	# Number of points for distribution
```"

# ╔═╡ 75955640-819a-11eb-3cee-57311a500d27
xcamthick = coordinates_to_camber_thickness(airfoil, 60)

# ╔═╡ 759667b0-819a-11eb-3301-07897e9ad255
md"You can do the inverse by feeding the abscissa, camber, and thickness distributions to:

```julia
camber_thickness_to_coordinates(xs, camber, thickness)
```
"

# ╔═╡ 7599c310-819a-11eb-2600-bd3ad40dd086
coords = camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3]);

# ╔═╡ 75a97a80-819a-11eb-18a7-cd074f39e623
md"You can split the coordinates into their upper and lower surfaces by calling:"

# ╔═╡ 75ab9d5e-819a-11eb-25f8-d5afd5b178e1
upper, lower = split_surface(airfoil);

# ╔═╡ 75b25420-819a-11eb-34d7-7b01e05a2dd2
md"### Plotting"

# ╔═╡ 75b895b0-819a-11eb-3094-458982d71f4b
x_upper, y_upper = getindex.(upper, 1), getindex.(upper, 2);

# ╔═╡ 75bd9ec0-819a-11eb-3eeb-0f0c7ad2ee27
x_lower, y_lower = getindex.(lower, 1), getindex.(lower, 2);

# ╔═╡ 75c4cab0-819a-11eb-3edc-f914f2239e9c
md"Let's plot everything!"

# ╔═╡ 755985c0-819a-11eb-1d29-b9dd9370799b
begin
	plot(aspectratio = 1)
	plot!(x_upper, y_upper, 
		  marker = :dot, markersize = 2, label = "NACA $(digits...) Upper")
	plot!(x_lower, y_lower, 
		  marker = :dot, markersize = 2, label = "NACA $(digits...) Lower")
	plot!(xcamthick[:,1], xcamthick[:,2], 
		  marker = :dot, markersize = 2, label = "NACA $(digits...) Camber")
	plot!(xcamthick[:,1], xcamthick[:,3], 
		  marker = :dot, markersize = 2, label = "NACA $(digits...) Thickness")
end

# ╔═╡ 2e9f7f40-819f-11eb-1a50-e7bb7cda3eab
md"""### Kulfan Class Shape Transformation Airfoils

AeroMDAO also provides the Class Shape Transformation method using Bernstein polynomials developed by Brenda Kulfan. You can read more about it here: [http://smtp.brendakulfan.com/docs/CST6.pdf](http://smtp.brendakulfan.com/docs/CST6.pdf)

In this parametrization, you define two arrays (of possible different lengths) which effectively determine the displacements at chord-wise locations, a tuple for the displacements of the upper and lower points at the trailing edge, an optional coefficient for leading edge modifications (provided in the reference above), and an optional number of points.

```julia
kulfan_CST(alpha_u      :: Vector{Real},    # Upper surface parameters
           alpha_l      :: Vector{Real},    # Lower surface parameters
           dzs          :: NTuple{2, Real}, # Upper and lower trailing edge points
           coeff_LE = 0 :: Real,            # Leading-edge modification coefficient
           n = 40       :: Integer)         # Number of points on each surface
```

"""

# ╔═╡ d607d6f0-81af-11eb-327f-d146fc83da3e
md"Try modifying the arrays below and see how the airfoil plot changes."

# ╔═╡ 3980af10-819f-11eb-09d8-779244f75746
begin
	alpha_u = [0.4, 0.3, 0.3, 0.15, 0.2]	
	alpha_l = [-0.4, -0.1, -0.07, -0.001]	
	dzs     = (0., 0.)						
	kulfan_airfoil = kulfan_CST(alpha_u, alpha_l, dzs, 0.2, 60)
end;

# ╔═╡ 6b2f1640-81af-11eb-0377-171002fee886
cst_upper, cst_lower = split_surface(kulfan_airfoil);

# ╔═╡ 75cc44be-819a-11eb-0815-098288453926
md"## Potential Flow Analysis

Now we have an airfoil, and we would like to analyze its aerodynamic characteristics. The potential flow panel method for inviscid analyses of airfoils, which you may have studied in your course on aerodynamics, provides decent estimations of the lift generated by the airfoil."

# ╔═╡ 75cd2f20-819a-11eb-329a-4539ebf53e84
md"To analyze our airfoil, we must convert the coordinates into a `Foil` type defined in `AeroMDAO`, as shown:"

# ╔═╡ 75d7b670-819a-11eb-0366-bf27460f6688
kulfan_foil = Foil(kulfan_airfoil, "Kulfan")

# ╔═╡ 75d8c7e0-819a-11eb-2eff-b512d6389fda
md"Our analysis also requires boundary conditions, which is the freestream flow defined by a magnitude ``V_\infty`` and angle of attack ``\alpha``. We provide these to the analysis by defining variables and feeding them to a `Uniform2D` type, corresponding to uniform flow in 2 dimensions."

# ╔═╡ 75e460a2-819a-11eb-1d01-27f5ab5a9e71
begin
	V = 1.0
	alpha = 5.0
end;

# ╔═╡ 75e59920-819a-11eb-1985-7180573a8376
uniform = Uniform2D(V, alpha)

# ╔═╡ 75f29170-819a-11eb-34df-551a908f9710
md"Now that we have our airfoil and boundary conditions, we can call the `solve_case()` function, which in this case has an associated method with the specification of ``n`` panels given by the optional argument `num_panels`, which is ``60`` by default. This will run the analysis and return the lift coefficient, the sectional lift and pressure coefficients, and the panels for the given case.
```julia
solve_case(foil 			:: Foil,
		   uniform 			:: Uniform2D; 
		   num_panels = 60 	:: Integer)
```
"

# ╔═╡ 092f3fe1-46fe-41c9-a4f0-f98294223bb1
naca_foil = Foil(airfoil, "NACA $digits")

# ╔═╡ 75fc5570-819a-11eb-1c77-4777dd857175
cl, cls, cms, cps, panels = solve_case(naca_foil, 
                                       uniform;
                                       viscous = false,
                                       sources = false, 
                                       wake_length = 1e5,
                                       wake_panels = 100,
                                       num_panels = 80);

# ╔═╡ 5fa7581e-81ad-11eb-2241-cd45f25bfec8
md"Note the difference between the lift coefficient computed and the sum of the sectional lift coefficients; this is due to numerical errors in the solution procedure."

# ╔═╡ 75fd18c0-819a-11eb-16c0-916dcd0b63a4
cl, sum(cls)

# ╔═╡ 760b70a0-819a-11eb-3888-291d1706ecde
md"### Post-processing
AeroMDAO provides more helper functions for post-processing data
"

# ╔═╡ 7f46182e-81a1-11eb-3dd3-6199e26a9ce2
begin
	tangents = panel_tangent.(panels)	# Tangents
	normals = panel_normal.(panels)		# Normals
	locs = panel_location.(panels)		# Upper or lower surface
end;

# ╔═╡ 87f7e4e0-81a1-11eb-0aff-d9781907f6ea
md"""
Let's see what the pressure and lift distribution curves look like over the airfoil. First, we define a function to get the $x$-locations and their corresponding values for the upper or lower surface, given arrays of panels and values.
"""

# ╔═╡ ad58b03e-fbc1-4f53-9747-d75f2f90eaba
# get_surface_values(panels, vals, surf = "upper") = partition(x -> (panel_location ∘ first)(x) == surf, (collect ∘ zip)(panels, vals), x -> ((first ∘ p1 ∘ first)(x), last(x)))

# ╔═╡ 6761aaa0-819b-11eb-37ea-1b5bba8a85ae
begin
	cp_lower, cp_upper = get_surface_values(panels, cps, "lower")
	cl_lower, cl_upper = get_surface_values(panels, cls, "lower")
end

# ╔═╡ c53a3170-81ae-11eb-28ed-89d8f798a6fc
begin
	plot(aspectratio = 1)
	plot!(getindex.(cst_upper, 1), getindex.(cst_upper, 2), 
		  marker = :dot, markersize = 2, label = "CST Upper")
	plot!(getindex.(cst_lower, 1), getindex.(cst_lower, 2), 
		  marker = :dot, markersize = 2, label = "CST Lower")
end

# ╔═╡ 9851e530-819b-11eb-13f8-794fffe0451c
begin
	plot(marker = 2, yflip = true, xlabel="(x/c)", ylabel = "Cp")
	plot!(cp_upper, label = "Upper", marker = :dot)
	plot!(cp_lower, label = "Lower", marker = :dot)
end

# ╔═╡ 676fdb70-819b-11eb-1d5b-d723542dab36
begin
	plot(marker = 2, xlabel="(x/c)", ylabel = "Cl")
	plot!(cl_upper, label = "Upper", marker = :dot)
	plot!(cl_lower, label = "Lower", marker = :dot)
end

# ╔═╡ 91025d10-819a-11eb-1664-33bca2c15831
begin
	alert(text) = Markdown.MD(Markdown.Admonition("warning", "Alert!", [text]))
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));
	definition(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("warning", "Definition", [text]));
	exercise(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("correct", "Exercise", [text]));
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
end;

# ╔═╡ 755f521e-819a-11eb-1e77-131f637a1454
alert(md"AeroMDAO uses the `StaticArrays` package to maximize computational speed. These share the same syntax as arrays, but you cannot mutate them (_i.e._ change their values).")

# ╔═╡ 760a8642-819a-11eb-3283-e749305d2c62
alert(md"""Support for drag prediction with boundary layer calculations will be added soon. For now, try out the amazing [Webfoil](http://webfoil.engin.umich.edu/) developed by the MDOLab at University of Michigan -- Ann Arbor!""")

# ╔═╡ 76325990-819a-11eb-170f-8144db240b6a
exercise(md"Run an analysis for a NACA-4 airfoil of your choice between angles ``-5`` and ``5`` with any increment you like and plot the lift coefficient variation against angle of attack.")

# ╔═╡ 763c92c0-819a-11eb-0e94-ad296f0a334a
hint(md"You can use ranges to generate the array of angles, and use a loop to execute the cases.")

# ╔═╡ Cell order:
# ╟─993a3ed2-819a-11eb-0464-ab84db58a8a2
# ╠═75407f80-819a-11eb-3307-e5c7d2dcef46
# ╠═7542a260-819a-11eb-36b0-4d6f6497e611
# ╠═75665700-819a-11eb-04e3-23d133e8dc6c
# ╠═526c4d80-f348-448a-a8f0-ff073bd96418
# ╟─754bca22-819a-11eb-0039-59eb17946634
# ╟─2060cc90-819f-11eb-1b1e-d1f61e2c878a
# ╟─754cdb90-819a-11eb-30f7-0d172a955548
# ╟─754dc5ee-819a-11eb-373d-f909ce6e2e3a
# ╠═755455a2-819a-11eb-1384-a184f1a5ed93
# ╠═7554f1e0-819a-11eb-06bb-6587881bb0a1
# ╟─755f521e-819a-11eb-1e77-131f637a1454
# ╟─758a32b0-819a-11eb-385d-e75b23f41ac4
# ╟─758b6b30-819a-11eb-35c6-99ad31016f75
# ╠═75955640-819a-11eb-3cee-57311a500d27
# ╟─759667b0-819a-11eb-3301-07897e9ad255
# ╠═7599c310-819a-11eb-2600-bd3ad40dd086
# ╟─75a97a80-819a-11eb-18a7-cd074f39e623
# ╠═75ab9d5e-819a-11eb-25f8-d5afd5b178e1
# ╟─75b25420-819a-11eb-34d7-7b01e05a2dd2
# ╠═75b895b0-819a-11eb-3094-458982d71f4b
# ╠═75bd9ec0-819a-11eb-3eeb-0f0c7ad2ee27
# ╟─75c4cab0-819a-11eb-3edc-f914f2239e9c
# ╠═755985c0-819a-11eb-1d29-b9dd9370799b
# ╟─2e9f7f40-819f-11eb-1a50-e7bb7cda3eab
# ╟─d607d6f0-81af-11eb-327f-d146fc83da3e
# ╠═3980af10-819f-11eb-09d8-779244f75746
# ╠═6b2f1640-81af-11eb-0377-171002fee886
# ╟─75cc44be-819a-11eb-0815-098288453926
# ╟─75cd2f20-819a-11eb-329a-4539ebf53e84
# ╠═75d7b670-819a-11eb-0366-bf27460f6688
# ╟─75d8c7e0-819a-11eb-2eff-b512d6389fda
# ╠═75e460a2-819a-11eb-1d01-27f5ab5a9e71
# ╠═75e59920-819a-11eb-1985-7180573a8376
# ╟─75f29170-819a-11eb-34df-551a908f9710
# ╠═092f3fe1-46fe-41c9-a4f0-f98294223bb1
# ╠═75fc5570-819a-11eb-1c77-4777dd857175
# ╟─5fa7581e-81ad-11eb-2241-cd45f25bfec8
# ╠═75fd18c0-819a-11eb-16c0-916dcd0b63a4
# ╟─760a8642-819a-11eb-3283-e749305d2c62
# ╟─760b70a0-819a-11eb-3888-291d1706ecde
# ╠═7f46182e-81a1-11eb-3dd3-6199e26a9ce2
# ╟─87f7e4e0-81a1-11eb-0aff-d9781907f6ea
# ╠═ad58b03e-fbc1-4f53-9747-d75f2f90eaba
# ╠═6761aaa0-819b-11eb-37ea-1b5bba8a85ae
# ╟─c53a3170-81ae-11eb-28ed-89d8f798a6fc
# ╟─9851e530-819b-11eb-13f8-794fffe0451c
# ╟─676fdb70-819b-11eb-1d5b-d723542dab36
# ╟─76325990-819a-11eb-170f-8144db240b6a
# ╟─763c92c0-819a-11eb-0e94-ad296f0a334a
# ╟─91025d10-819a-11eb-1664-33bca2c15831
# ╟─91040ac0-819a-11eb-2967-2909f8487f5b
