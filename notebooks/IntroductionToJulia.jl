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

<center><img src=https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/b5551ca7946b4a25746c045c15fbb8806610f8d0/images/julia-logo-color.svg></center>"""

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

```bash
julia> using Pluto
julia> Pluto.run()
```
"""

# ╔═╡ f0111d90-5f09-11eb-3977-237d140094d2
md"## First steps"

# ╔═╡ 24808810-5eef-11eb-0d17-812e03dcd966
md"""
### Running your first program

```bash
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

# ╔═╡ 7cc1eeb0-5f0d-11eb-39c6-2715c40576b0
xs = [1, 2, 3, 4.]

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
md"1. Expressions"

# ╔═╡ 4faa88e0-5f0a-11eb-1d42-c51848513150
begin
	j = 0
	for i in 1:10
		j += i
	end
	j
end

# ╔═╡ 9408abc0-5f0a-11eb-35f2-cd1b2fba577b
md"2. Comprehensions (usually nicer and more readable)"

# ╔═╡ 9f4d3140-5f0a-11eb-3302-8f26d4b0a6b7
m = sum(i for i in 1:10)

# ╔═╡ a0ed1080-5f08-11eb-0677-215e424d8299
md"### Your first function"

# ╔═╡ 7b4f39f0-5f0f-11eb-097d-cff9078fe957
md"""
Let's try to generalise our previous computations for generic $n$ natural numbers, by using _functions_ which take arguments and return values.

There are three ways to write functions:
""" 

# ╔═╡ 1b524420-5f0a-11eb-3613-fd56bcd66dd5


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
	correct(text=md"Great! You got the right answer! Let's move on to the next section.") = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]));
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]));
end;

# ╔═╡ 108e4700-5f0f-11eb-34f3-2509225fb880
hint(md"This will still work if you remove the braces, and it'll be more efficient! This usage is called a _generator_ expression.")

# ╔═╡ Cell order:
# ╟─d0e4dac0-5edb-11eb-0312-678a24ad0573
# ╟─0f4bdbc0-5eea-11eb-35c3-d13c5d187b11
# ╟─5f6f01f0-5eee-11eb-090f-07ee55921d29
# ╟─dc211a80-5ee9-11eb-2ec4-25e7f37d0183
# ╟─43c37610-5eef-11eb-06ed-fdc3053b6c20
# ╠═f0111d90-5f09-11eb-3977-237d140094d2
# ╟─24808810-5eef-11eb-0d17-812e03dcd966
# ╟─1e041650-5f0d-11eb-2c1d-69a91cd398ae
# ╠═2509e0b0-5f0d-11eb-0b86-ef2bfdc59221
# ╠═76494920-5f0d-11eb-040e-63b283f2f031
# ╠═7cc1eeb0-5f0d-11eb-39c6-2715c40576b0
# ╟─ed0d2c60-5f09-11eb-3346-73bcb62dc713
# ╠═3e238902-5f0a-11eb-27e9-61eb774d3620
# ╠═2309d252-5f0a-11eb-3db5-d16147bb8329
# ╟─e56e2f40-5f09-11eb-29a9-2db4ee699115
# ╟─591a2752-5f0a-11eb-3110-d17376a0e3a9
# ╟─6d741d50-5f0a-11eb-3da2-e977918b03c2
# ╠═4faa88e0-5f0a-11eb-1d42-c51848513150
# ╟─9408abc0-5f0a-11eb-35f2-cd1b2fba577b
# ╠═9f4d3140-5f0a-11eb-3302-8f26d4b0a6b7
# ╟─108e4700-5f0f-11eb-34f3-2509225fb880
# ╟─a0ed1080-5f08-11eb-0677-215e424d8299
# ╟─7b4f39f0-5f0f-11eb-097d-cff9078fe957
# ╠═1b524420-5f0a-11eb-3613-fd56bcd66dd5
# ╟─efa35d10-5f08-11eb-2108-8373143309cb
# ╟─ee5d60d0-5f0e-11eb-0574-6553f2bad5a5
