### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ c24b1700-9134-11eb-2261-75af42b46ddc
begin
	using LinearAlgebra
	using Plots
	using ControlSystems
	using DifferentialEquations
end

# ╔═╡ ae449830-9134-11eb-13ca-8fe9c4eb6544
md"# Topics in Nonlinear Dynamics"

# ╔═╡ cc46efe0-9134-11eb-05c7-8f440337503d
md"""
Dynamical systems are defined as a set of quantities, possibly interdependent, which are all dependent on time. This is written as a 'non-autonomous' system, with $\mathbf f\colon \mathbb R^n \times \mathbb R \to \mathbb R^n$.

```math
\dot{\mathbf{x}} = \mathbf f(\mathbf x, t)
```

Technically, you can define $t = x_{n+1}$ and convert it into an autonomous system, with supposedly difficult consequences in interpretation.
"""

# ╔═╡ 7e426a50-9138-11eb-0501-034f297f100a
begin
	f₁(x, t) = x[1]^2 - x[2]^3 + x[3] * t
	f₂(x, t) = x[1] * t + x[2] - x[3] * t^2
	f₃(x, t) = x[1] * x[2] + x[3] * t
	
	function f!(du, u, p, t)
		 du[1] = p[1]*(u[2]-u[1])
		 du[2] = u[1]*(p[2]-u[3]) - u[2]
		 du[3] = u[1]*u[2] - p[3]*u[3]
	end
end

# ╔═╡ f4a849e2-9146-11eb-2be4-0da2d510abb1
begin
	u0 = [1.0, 0.1, 0.1]
	p0 = [10, 28, 8/3]
	tspan = (0.0, 100.0)
	prob = ODEProblem(f!,u0,tspan,p0)
	sol = solve(prob)
end;

# ╔═╡ 6785a750-9147-11eb-0d7c-f754111c9ba8
plot(sol, vars=(1,2,3))

# ╔═╡ Cell order:
# ╟─ae449830-9134-11eb-13ca-8fe9c4eb6544
# ╠═c24b1700-9134-11eb-2261-75af42b46ddc
# ╟─cc46efe0-9134-11eb-05c7-8f440337503d
# ╠═7e426a50-9138-11eb-0501-034f297f100a
# ╠═f4a849e2-9146-11eb-2be4-0da2d510abb1
# ╠═6785a750-9147-11eb-0d7c-f754111c9ba8
