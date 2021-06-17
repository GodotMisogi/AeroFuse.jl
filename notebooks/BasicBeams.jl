### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ fa875bc6-f285-42fa-9166-f99ef04cb8c6
using AeroMDAO

# ╔═╡ 47729952-4ff5-4a7a-b3bd-6aa5be1a1e5e
md"### Stiffness matrix"

# ╔═╡ cd8447c0-c8e0-11eb-0430-fbea74d11fb5
# Deflection stiffness matrix
K = deflection_stiffness_matrix([1., 1.], [1., 1.], [2., 2.], :z)

# ╔═╡ 6bb731d4-5aed-4b1a-8c2c-60e2ccb94b77
md"### 1. Fixed hinged beam subjected to force and moment at the center"

# ╔═╡ b7b7a275-7ee3-42af-a351-f670dba51220
A1 = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3 

# ╔═╡ 61ba9229-ce69-4510-a3ec-099db04d075e
b1 = [-1000, 1000, 0]    # F2, M2, M3

# ╔═╡ 2c60a4b5-10f9-4902-9ded-cf49f3b46719
x1 = A1 \ b1

# ╔═╡ eb968bf4-868e-4935-afbc-c689e470c405
md"Forces"

# ╔═╡ 1809373e-26e3-4e97-98a8-2bbbe44a730a

F1 = K * [ 0.; 0.; x1[1:2]; 0.; x1[3] ]

# ╔═╡ b75d1bc9-362e-4432-88b2-2355e54fc581
md"### 2. Propped cantilever beam with force at one end"

# ╔═╡ 0947e16c-817a-45f6-a32a-6a78c2368861

A2 = K[[1,2,4],[1,2,4]] # v1, φ1, φ2

# ╔═╡ 0d1fed4b-ffb7-41c2-bde4-fe118a989fb5
b2 = [10, 0, 0]

# ╔═╡ 9a4305d8-2bd9-4c6f-99ca-4f6b82a29601
x2 = A2 \ b2

# ╔═╡ 1daccfa1-5553-4f6b-8555-0005a6d63dfa
## Forces
F2 = K * [ x2[1:2]; 0.; x2[3]; 0.; 0. ]

# ╔═╡ Cell order:
# ╠═fa875bc6-f285-42fa-9166-f99ef04cb8c6
# ╟─47729952-4ff5-4a7a-b3bd-6aa5be1a1e5e
# ╠═cd8447c0-c8e0-11eb-0430-fbea74d11fb5
# ╟─6bb731d4-5aed-4b1a-8c2c-60e2ccb94b77
# ╠═b7b7a275-7ee3-42af-a351-f670dba51220
# ╠═61ba9229-ce69-4510-a3ec-099db04d075e
# ╠═2c60a4b5-10f9-4902-9ded-cf49f3b46719
# ╟─eb968bf4-868e-4935-afbc-c689e470c405
# ╠═1809373e-26e3-4e97-98a8-2bbbe44a730a
# ╟─b75d1bc9-362e-4432-88b2-2355e54fc581
# ╠═0947e16c-817a-45f6-a32a-6a78c2368861
# ╠═0d1fed4b-ffb7-41c2-bde4-fe118a989fb5
# ╠═9a4305d8-2bd9-4c6f-99ca-4f6b82a29601
# ╠═1daccfa1-5553-4f6b-8555-0005a6d63dfa
