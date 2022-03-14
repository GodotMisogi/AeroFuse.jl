# AeroMDAO -- Aircraft Design Platform

AeroMDAO is meant to be a toolbox for aircraft design analyses. It currently provides convenient methods for developing studies in aerodynamics and structures, with implementations in other relevant fields such as flight dynamics and propulsion in progress.

**Authors**: Arjit Seth, Stephane Redonnet, and Rhea P. Liem

## Objectives

The current focus is to enable tutorials in computation for undergraduates in an aerospace educational curriculum, particularly at The Hong Kong University of Science and Technology. For this purpose, the code is written in a functional style replicating the mathematics presented in textbooks as much as possible.

An additional aim is to write code compatible with automatic differentiation libraries written in Julia.

## Features

**Disclaimer**: The implementations are work-in-progress, and hence the results may not be entirely accurate. Please exercise caution when interpreting the results until validation cases are added.

**Functioning:**

- Basic geometric tools for airfoil processing and wing design. 
- Steady, inviscid, incompressible, isentropic 2D and 3D aerodynamic analyses.
  - 2D - Doublet-source panel method
  - 3D - Vortex lattice method with derivatives
- Empirical viscous drag prediction methods for wings.
- Finite-element beam model based on Euler-Bernoulli theory for structural analysis.

![](assets/VortexLattice.svg)
![](assets/LinearVortex.svg)

**In progress:**

- Aeroelastic analyses coupling vortex lattice and beam element models.
- 2D integral boundary layer solver for viscous-inviscid analyses.
- Generic, flexible and efficient flight dynamics integrator in 2 and 3 dimensions.
- Fuselage geometry parametrization and aerodynamic implementation.

![](https://godot-bloggy.xyz/post/diagrams/AerostructAircraft.svg)

## Installation

```julia
julia> using Pkg; Pkg.add("AeroMDAO")
julia> Pkg.test("AeroMDAO")
julia> using AeroMDAO
```

## Citation

If you use AeroMDAO in your research, please cite the following until any relevant material is actually published:

```bibtex
@software{aeromdao,
  author  = {Arjit Seth, Rhea P. Liem, Stephane Redonnet},
  title   = {AeroMDAO},
  url     = {https://github.com/GodotMisogi/AeroMDAO},
  version = {0.3.7},
  date    = {2021-03-14},
}
```

## Contents

```@contents
Pages = ["geometry.md", "aerodynamics.md", "structures.md", "theory.md", "api.md"]
Depth =  5
```