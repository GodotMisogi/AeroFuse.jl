### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ d0e4dac0-5edb-11eb-0312-678a24ad0573
md"""
# MECH 3620 -- Introduction to Julia
"""

# ╔═╡ 0f4bdbc0-5eea-11eb-35c3-d13c5d187b11
html"""
<center><img src=https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/b5551ca7946b4a25746c045c15fbb8806610f8d0/images/julia-logo-color.svg></center>
"""

# ╔═╡ 5f6f01f0-5eee-11eb-090f-07ee55921d29
md"""## Why Julia?
> Most (insert branch of engineering here) students and instructors are familiar with MATLAB.
1. Its syntax is very similar to MATLAB's.
2. It's faster than MATLAB.
3. It has Pluto, this nice notebook format for running experiments and giving instructive presentations.
"""

# ╔═╡ dc211a80-5ee9-11eb-2ec4-25e7f37d0183
md"""
## Installation

Go to [https://julialang.org/downloads](https://julialang.org/downloads) and download the current stable release, using the correct version for your operating system (Linux x86, Mac, Windows, etc).

---

Cheatsheets: 

[MATLAB-Julia-Python comparative cheatsheet by QuantEcon group
](https://cheatsheets.quantecon.org/)

[Fastrack to Julia](https://juliadocs.github.io/Julia-Cheat-Sheet/)
"""

# ╔═╡ 43c37610-5eef-11eb-06ed-fdc3053b6c20
md"""
## Running a Pluto notebook

```julia
julia> using Pluto
julia> Pluto.run()
```
"""

# ╔═╡ f0111d90-5f09-11eb-3977-237d140094d2
md"## First steps"

# ╔═╡ 24808810-5eef-11eb-0d17-812e03dcd966
md"""
### Running your first program

```julia
julia> println("Hello world!")
Hello world!
```
"""

# ╔═╡ 1e041650-5f0d-11eb-2c1d-69a91cd398ae
md"### Your first variable"

# ╔═╡ 2509e0b0-5f0d-11eb-0b86-ef2bfdc59221
david = 0.5

# ╔═╡ 76494920-5f0d-11eb-040e-63b283f2f031
md"### Your first array"

# ╔═╡ a3f4b6f0-5f23-11eb-09ba-0bed5674f6a9
md"1D array:"

# ╔═╡ 7cc1eeb0-5f0d-11eb-39c6-2715c40576b0
xs = [1, 2, 3, 4]

# ╔═╡ a88d2ad0-5f23-11eb-2e9b-45a2858e2db3
md"2D array:"

# ╔═╡ ac1c9320-5f23-11eb-1e76-cdcbc766ecec
ys = [ 1 2 3 ; 
	   4 5 6.]

# ╔═╡ 467059e0-5f2c-11eb-3d82-ed166a16040c
md"You can also transpose arrays:"

# ╔═╡ 50c32c60-5f2c-11eb-1130-2f9ceb28caea
ys'

# ╔═╡ 09aae9c0-5fb4-11eb-2897-2dabe0b66efd
md"Notice the different types mentioned above, e.g. `LinearAlgebra.Adjoint`. We will discuss type systems shortly."

# ╔═╡ ed0d2c60-5f09-11eb-3346-73bcb62dc713
md"### Your first conditional"

# ╔═╡ 3e238902-5f0a-11eb-27e9-61eb774d3620
cuck = false

# ╔═╡ 2309d252-5f0a-11eb-3db5-d16147bb8329
if !cuck
	return johnny(denny, mark)
else
	return lisa(mark, claudette)
end

# ╔═╡ e56e2f40-5f09-11eb-29a9-2db4ee699115
md"### Your first loop"

# ╔═╡ 591a2752-5f0a-11eb-3110-d17376a0e3a9
md"There are two ways to write `for` loops:"

# ╔═╡ 6d741d50-5f0a-11eb-3da2-e977918b03c2
md"1. Expressions:"

# ╔═╡ 4faa88e0-5f0a-11eb-1d42-c51848513150
begin
	j = 0
	for i in 1:10
		j += i
	end
	j
end

# ╔═╡ 9408abc0-5f0a-11eb-35f2-cd1b2fba577b
md"2. Comprehensions (usually nicer and more readable):"

# ╔═╡ 4eb94c3e-5f24-11eb-0479-a75fa79801d2
md"Note that this operation of summation is a common operation, and you could use Julia's `sum` function which generically sums the elements of an array of numbers."

# ╔═╡ 9f4d3140-5f0a-11eb-3302-8f26d4b0a6b7
m = sum([ i for i in 1:10 ])

# ╔═╡ a0ed1080-5f08-11eb-0677-215e424d8299
md"### Your first function"

# ╔═╡ 7b4f39f0-5f0f-11eb-097d-cff9078fe957
md"""
Let's try to generalise our previous computations for generic $n$ natural numbers, by using _functions_ which take arguments and return values.

There are three ways to write functions:
""" 

# ╔═╡ 1b524420-5f0a-11eb-3613-fd56bcd66dd5
md"1. 'Computer Science' syntax:"

# ╔═╡ cb6cc78e-5f23-11eb-25a1-5db7e3643ceb
function sum_nat_num(n :: Integer)
	return sum(1:n)   
end

# ╔═╡ 21322d50-5f24-11eb-3e8e-092d8a5d1dfc
md"The `return` keyword is optional; Julia function blocks return the value of their last statement."

# ╔═╡ 3bd24730-5f24-11eb-34b6-b92b436ff3f0
sum_nat_num(10)

# ╔═╡ d9512540-5f23-11eb-32f6-eb121a922fb7
md"2. 'Mathematical' syntax (preferred for readability):"

# ╔═╡ e0f1d010-5f23-11eb-2146-5fd2ca7a3e40
sum_nat_num_math(n :: Integer) = sum(1:n)

# ╔═╡ 40067b00-5f24-11eb-1dec-611fc8bd6736
sum_nat_num_math(10)

# ╔═╡ 4f737f0e-5f25-11eb-305e-fd2c5053f59a
md"3. ``\lambda``-expressions:"

# ╔═╡ 5d83f5d0-5f25-11eb-0a9b-890e2b3f4ff2
sum_nat_num_λ = n :: Integer -> sum(1:n)

# ╔═╡ 661468c0-5f2e-11eb-1272-774164f1904d
md"Notice how `sum_nat_num_λ` is now denoted as a _variable_ with a function assigned to it as a value. However, you cannot change its value. You can still use it like a function."

# ╔═╡ 5903aec0-5f2e-11eb-2411-25a29958ce46
sum_nat_num_λ(10)

# ╔═╡ efa35d10-5f08-11eb-2108-8373143309cb
html"""<style>
main {
    max-width: 60%;
}
"""

# ╔═╡ ee5d60d0-5f0e-11eb-0574-6553f2bad5a5
begin
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]));
	alert(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Alert!", [text]));
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
end;

# ╔═╡ 108e4700-5f0f-11eb-34f3-2509225fb880
hint(md"This will still work if you remove the braces, and it'll be more efficient! 	
```julia
	m = sum(i for i in 1:10)
```
This usage is called a _generator_ expression. Copy it to a new cell and try it out!")

# ╔═╡ 94f0f4a0-5f25-11eb-20b1-3d1cf2900c0b
alert(
md"""``\lambda``-expressions have a specific usage, and are not encouraged for writing generic functions. 
	
This specific usage is for defining temporary functions to perform small operations that do not require functions of their own generically.

e.g. Summing an array of squares of natural numbers:
```julia
	sum(x -> x^2, 1:10)
```""")

# ╔═╡ Cell order:
# ╟─d0e4dac0-5edb-11eb-0312-678a24ad0573
# ╟─0f4bdbc0-5eea-11eb-35c3-d13c5d187b11
# ╟─5f6f01f0-5eee-11eb-090f-07ee55921d29
# ╟─dc211a80-5ee9-11eb-2ec4-25e7f37d0183
# ╟─43c37610-5eef-11eb-06ed-fdc3053b6c20
# ╟─f0111d90-5f09-11eb-3977-237d140094d2
# ╟─24808810-5eef-11eb-0d17-812e03dcd966
# ╟─1e041650-5f0d-11eb-2c1d-69a91cd398ae
# ╠═2509e0b0-5f0d-11eb-0b86-ef2bfdc59221
# ╟─76494920-5f0d-11eb-040e-63b283f2f031
# ╟─a3f4b6f0-5f23-11eb-09ba-0bed5674f6a9
# ╠═7cc1eeb0-5f0d-11eb-39c6-2715c40576b0
# ╟─a88d2ad0-5f23-11eb-2e9b-45a2858e2db3
# ╠═ac1c9320-5f23-11eb-1e76-cdcbc766ecec
# ╠═467059e0-5f2c-11eb-3d82-ed166a16040c
# ╠═50c32c60-5f2c-11eb-1130-2f9ceb28caea
# ╟─09aae9c0-5fb4-11eb-2897-2dabe0b66efd
# ╟─ed0d2c60-5f09-11eb-3346-73bcb62dc713
# ╠═3e238902-5f0a-11eb-27e9-61eb774d3620
# ╠═2309d252-5f0a-11eb-3db5-d16147bb8329
# ╟─e56e2f40-5f09-11eb-29a9-2db4ee699115
# ╟─591a2752-5f0a-11eb-3110-d17376a0e3a9
# ╟─6d741d50-5f0a-11eb-3da2-e977918b03c2
# ╠═4faa88e0-5f0a-11eb-1d42-c51848513150
# ╟─9408abc0-5f0a-11eb-35f2-cd1b2fba577b
# ╟─4eb94c3e-5f24-11eb-0479-a75fa79801d2
# ╠═9f4d3140-5f0a-11eb-3302-8f26d4b0a6b7
# ╟─108e4700-5f0f-11eb-34f3-2509225fb880
# ╟─a0ed1080-5f08-11eb-0677-215e424d8299
# ╟─7b4f39f0-5f0f-11eb-097d-cff9078fe957
# ╟─1b524420-5f0a-11eb-3613-fd56bcd66dd5
# ╠═cb6cc78e-5f23-11eb-25a1-5db7e3643ceb
# ╟─21322d50-5f24-11eb-3e8e-092d8a5d1dfc
# ╠═3bd24730-5f24-11eb-34b6-b92b436ff3f0
# ╟─d9512540-5f23-11eb-32f6-eb121a922fb7
# ╠═e0f1d010-5f23-11eb-2146-5fd2ca7a3e40
# ╠═40067b00-5f24-11eb-1dec-611fc8bd6736
# ╟─4f737f0e-5f25-11eb-305e-fd2c5053f59a
# ╠═5d83f5d0-5f25-11eb-0a9b-890e2b3f4ff2
# ╟─661468c0-5f2e-11eb-1272-774164f1904d
# ╠═5903aec0-5f2e-11eb-2411-25a29958ce46
# ╟─94f0f4a0-5f25-11eb-20b1-3d1cf2900c0b
# ╟─efa35d10-5f08-11eb-2108-8373143309cb
# ╟─ee5d60d0-5f0e-11eb-0574-6553f2bad5a5
