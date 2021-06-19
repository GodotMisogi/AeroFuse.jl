# AeroMDAO

AeroMDAO is meant to be a toolbox for aircraft design analyses. It currently provides convenient methods for developing studies in aerodynamics, with aims to develop implementations in other relevant fields such as structures, propulsion, stability, etc.

## Aims

The current focus is to enable tutorials in computation for undergraduates in an aerospace educational curriculum, particularly at The Hong Kong University of Science and Technology. For this purpose, the code is written in a functional style replicating the mathematics presented in textbooks as much as possible.

**Aside**: This also turns out to develop very efficient (and more readable!) code without the need for functions that modify their arguments in-place (read mutation), which is against functional style. Of course, the code is not written in a higher-order Haskellian style using monads or any other category theory constructs, although that would probably be very interesting and pleasing to implement! *Hint*: Lists are secretly monads.

## Features

AeroMDAO currently provides basic geometric tools for airfoil processing, and panel methods for inviscid 2D analyses. A vortex lattice method presented in Mark Drela's *Flight Vehicle Aerodynamics* has also been implemented for inviscid 3D analyses with some traditional viscous drag prediction methods.

**Disclaimer**: The current implementation is a major work-in-progress, and hence the results may not be entirely accurate. It has extensively avoided referring to other implementations for originality. Please exercise caution when interpreting the results until validation cases are added.

## Installation

```julia
julia> using Pkg; Pkg.add("AeroMDAO")
julia> Pkg.test("AeroMDAO")
julia> using AeroMDAO
```

## Examples

Refer to the [cases document](examples/cases.md) in the `examples` folder to see how to get started with aerodynamic analyses using AeroMDAO.

## Citation

If you use AeroMDAO in your research, please cite the following until any relevant material is actually published:

```bibtex
@software{aeromdao,
  author  = {Arjit Seth, Rhea P. Liem, Stephane Redonnet},
  title   = {AeroMDAO},
  url     = {https://github.com/GodotMisogi/AeroMDAO},
  version = {0.2.3},
  date    = {2021-06-01},
}
```