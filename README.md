# AeroMDAO

AeroMDAO is meant to be a toolbox for aircraft design analyses. It currently provides convenient methods for developing studies in aerodynamics, with aims to develop implementations in other relevant fields such as structures, propulsion, stability, etc.

## Installation

```julia
julia> using Pkg; Pkg.add("https://github.com/GodotMisogi/AeroMDAO")
julia> Pkg.test(AeroMDAO)
julia> using AeroMDAO
```
---

## Aircraft Geometry

### Airfoil Parametrization

*NACA 4-digit Airfoil Parametrization*

```julia
airfoil = naca4(digits :: NTuple{4, Real})
```

```julia
airfoil = naca4((2,4,1,2))
```

*Kulfan Class Shape Transformation (CST) Method*
```julia
foil = kulfan_CST(alpha_u       :: Vector{Real}, 
                  alpha_l       :: Vector{Real},
                  dzs           :: NTuple{2, Real},
                  coeff_LE = 0. :: Real,
                  n = 0         :: Integer)
```

```julia
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
dzs = (1e-4, 1e-4)
foil = kulfan_CST(alpha_u, alpha_l, dzs, 0.2)
```

*Camber-Thickness Representation*
```julia
foil_camthick(coords    :: Array{2, Real}, 
              num = 40  :: Integer)
```

```julia
foilpath = "path/to/your/airfoil.dat"   # Airfoil coordinates file path
coords = read_foil(foilpath)            # Read coordinates file
cos_foil = cosine_foil(coords, 51)      # Cosine spacing with 51 points on upper and lower surfaces
xcamthick = foil_camthick(cos_foil)     # Convert to camber-thickness representation
foiler = camthick_foil(xcamthick[:,1],  # x-components
                       xcamthick[:,2],  # Camber distribution
                       xcamthick[:,3])  # Thickness distribution
```

### Doublet-Source Panel Method

AeroMDAO provides convenience functions using its specific types for analyses. To analyse this airfoil using the doublet-source panel method, the following method is called:

```julia
solve_case(foil             :: Foil, 
           uniform          :: Uniform2D;
           num_panels = 60  :: Integer)
```

```@example aeromdao
airfoil = naca4((2,4,1,2)) # hide
V, α = 1.0, 3.0 

foil = Foil(airfoil)
uniform = Uniform2D(V, α)
cl = solve_case(foil, uniform, num_panels = 60)
```

### Wing Parametrization

The following image depicts the parametrization schema used for wing planforms in terms of nested trapezoids.

![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

This parametrization is implemented as a composite type called a `HalfWing`: 

```julia
HalfWing(foils  	:: Vector{Foil}, 	# Foil profiles
         chords 	:: Vector{Real}, 	# Chord lengths
         twists 	:: Vector{Real}, 	# Twist angles (deg)
         spans  	:: Vector{Real}, 	# Section span lengths
         dihedrals	:: Vector{Real},	# Dihedral angles (deg)
         sweeps 	:: Vector{Real})	# Leading-edge sweep angles (deg)
```

We can create a `Wing` by feeding two `HalfWing`s to it:
```julia
Wing(left :: HalfWing, right :: HalfWing)
```

_e.g._
```julia
airfoil = naca4((2,4,1,2))
foils = Foil.(airfoil for i in 1:3)
wing_right = HalfWing(foils,
                      [0.4, 0.2, 0.1],
                      [0., 2., 5.],
                      [1.0, 0.1],
                      [0., 60.],
                      [0., 30.])

wing = Wing(wing_right, wing_right)
```

### Vortex Lattice Method

For a 3D case, we use the vortex lattice method for initial designs, given its quick speed for fast analyses. For an analysis, you require ambient reference conditions. In this case, you need the density $\rho$ and a reference location $x_\text{ref}$ for calculating forces and moments.


Now we run the case with specifications of the number of spanwise and chordwise panels by calling the `solve_case()` function, which has an associated method:
```julia
solve_case(wing                     :: Union{Wing, HalfWing},
           freestream               :: Freestream, 
           ρ                        :: Real, 
           r_ref = [0.25, 0., 0.]; 
           span_num = 5             :: Integer, 
           chord_num = 10           :: Integer)
```

```julia
ρ = 1.225
r_ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
U = 10.0
α = 5.0
β = 0.0
Ω = [0.0, 0.0, 0.0]
freestream = Freestream(U, α, β, Ω);
nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, freestream, ρ, r_ref; span_num = 5, chord_num = 10);
```