# AeroMDAO

A toolbox for aircraft design analyses. (WIP)

## Installation

```julia
julia> using Pkg; Pkg.add("https://github.com/GodotMisogi/AeroMDAO")
julia> Pkg.test(AeroMDAO)
julia> using AeroMDAO
```
---

## Aircraft Geometry

### Your First Airfoil Analysis

*NACA 4-digit Airfoil Parametrization*:

```julia
airfoil = naca4(digits :: NTuple{4, Real})
```

*Doublet-Source Panel Method*:

AeroMDAO provides convenience functions using its specific types for analyses. To analyse this airfoil using the doublet-source panel method, the following method is called:

```julia
solve_case(foil :: Foil, uniform :: Uniform2D; num_panels :: Integer)
```
_e.g._
```julia
airfoil = naca4((2,4,1,2))
V, α = 1.0, 3.0 

foil = Foil(airfoil)
uniform = Uniform2D(V, α)
cl = solve_case(foil, uniform, num_panels = 60)
```

### Your First Wing

The following image depicts the parametrization schema used for wing planforms in terms of nested trapezoids.

![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

This parametrization is implemented as a composite type called a `HalfWing`: 

```julia
HalfWing(
	foils  	:: Vector{Foil}, 	# Foil profiles
	chords 	:: Vector{Real}, 	# Chord lengths
	twists 	:: Vector{Real}, 	# Twist angles (deg)
	spans  	:: Vector{Real}, 	# Section span lengths
	dihedrals  :: Vector{Real},	# Dihedral angles (deg)
	sweeps 	:: Vector{Real}		# Leading-edge sweep angles (deg)
	)
```

We can create a `Wing` by feeding two `HalfWing`s to it:
```julia
Wing(left :: HalfWing, right :: HalfWing)
```

_e.g._
```julia
airfoil = naca4((2,4,1,2))
foils = Foil.(airfoil for i in 1:3)
wing_right = HalfWing(
    foils,
    [0.4, 0.2, 0.1],
    [0., 2., 5.],
    [1.0, 0.1],
    [0., 60.],
    [0., 30.]
    )

wing = Wing(wing_right, wing_right)
```