
```@setup aeromdao
using AeroMDAO
using Plots
using LaTeXStrings
gr()
```

The geometry setup for aircraft analyses is described here.

## Airfoil Parametrization

AeroMDAO provides some basic parametrizations commonly used for airfoils, with conversion to the camber-thickness representation and vice versa.

*NACA 4-digit Airfoil Parametrization*

```julia
naca4(digits :: NTuple{4, Real};  # Digits, e.g. (2,4,1,2)
      sharp_trailing_edge = true) # Sharp or blunt trailing edge
```

```julia
airfoil = naca4(2,4,1,2)
```

*Kulfan Class Shape Transformation (CST) Method*
```julia
kulfan_CST(alpha_u                 :: Vector{Real},    # Upper surface parameters
           alpha_l                 :: Vector{Real},    # Lower surface parameters
           (dz_u, dz_l) = (0., 0.) :: NTuple{2, Real}, # Upper and lower trailing edge points
           (LE_u, LE_l) = (0., 0.) :: NTuple{2, Real}, # Leading-edge modification coefficient
           n = 40                  :: Integer)         # Number of points on each surface
```

```julia
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
dzs     = (1e-4, 1e-4)
foil    = kulfan_CST(alpha_u, alpha_l, dzs, 0.2)
```

```@example aeromdao
digits  = (2,4,1,2)
airfoil = naca4((digits))
xcamthick = camber_thickness(airfoil, 60)
coords = camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])
upper, lower = split_surface(airfoil);
x_upper, y_upper = upper[:,1], upper[:,2]
x_lower, y_lower = lower[:,1], lower[:,2]
plot(aspectratio = 1)
af_plot = plot(aspect_ratio = 1, xlabel = L"(x/c)", ylabel = L"y")
plot!(x_upper, y_upper, ls = :solid, lw = 2, c = :cornflowerblue, label = "NACA $(digits...) Upper")
plot!(x_lower, y_lower, ls = :solid, lw = 2, c = :orange, label = "NACA $(digits...) Lower")
plot!(xcamthick[:,1], xcamthick[:,2], ls = :dash, lw = 2, c = :burlywood3, label = "NACA $(digits...) Camber")
plot!(xcamthick[:,1], xcamthick[:,3], ls = :dash, lw = 2, c = :grey, label = "NACA $(digits...) Thickness")
```

## Wing Parametrization

### HalfWing

The following image depicts the parametrization schema used for wing planforms in terms of nested trapezoids.

![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

This parametrization is implemented as a composite type called a `HalfWing`.

```julia
HalfWing(foils     :: Vector{Foil}, # Foil profiles
         chords    :: Vector{Real}, # Chord lengths
         twists    :: Vector{Real}, # Twist angles (deg)
         spans     :: Vector{Real}, # Section span lengths
         dihedrals :: Vector{Real}, # Dihedral angles (deg)
         LE_sweeps :: Vector{Real}; # Leading-edge sweep angles (deg)
         position = zeros(3),
         angle    = 0.,
         axis     = [1., 0., 0.])
```

### Wing

We can create a `Wing` by feeding two `HalfWing`s to it. Alternatively, a convenience method is provided by simply replacing `HalfWing` with `Wing` using the previous definition.

```julia
Wing(foils     :: Vector{Foil}, # Foil profiles
     chords    :: Vector{Real}, # Chord lengths
     twists    :: Vector{Real}, # Twist angles (deg)
     spans     :: Vector{Real}, # Section span lengths
     dihedrals :: Vector{Real}, # Dihedral angles (deg)
     LE_sweeps :: Vector{Real}; # Leading-edge sweep angles (deg)
     position = zeros(3),
     angle    = 0.,
     axis     = [1., 0., 0.])
```

There is also a convenient function for printing this information, whose last optional argument provides the name of the wing. [Pretty Tables](https://github.com/ronisbr/PrettyTables.jl) is used for pretty-printing.

We can obtain the relevant geometric information of the wing for design analyses by calling convenience methods, which automatically perform the necessary calculations on the nested trapezoidal planforms.

```julia
AR = aspect_ratio(wing)
b  = span(wing)
S  = projected_area(wing)
c  = mean_aerodynamic_chord(wing)

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
```

#### Additional Methods

Wings with single sections can be defined with the following constructors and keywords.

```julia
HalfWingSection(; span, dihedral, LE_sweep, taper, root_chord,
                  root_twist, tip_twist, root_foil, tip_foil,
                  position, angle, axis)

WingSection(; span, dihedral, LE_sweep, taper, root_chord,
              root_twist, tip_twist, root_foil, tip_foil,
              position, angle, axis)
```
**Arguments**
- `span       :: Real         = 1.`: Span length 
- `dihedral   :: Real         = 1.`: Dihedral angle (degrees)
- `LE_sweep   :: Real         = 0.`: Leading-edge sweep angle (degrees)
- `taper      :: Real         = 1.`: Taper ratio of tip to root chord
- `root_chord :: Real         = 1.`: Root chord length
- `root_twist :: Real         = 0.`: Twist angle at root (degrees)
- `tip_twist  :: Real         = 0.`: Twist angle at tip (degrees)
- `root_foil  :: Array{Real}  = naca4((0,0,1,2))`: Foil coordinates at root
- `tip_foil   :: Array{Real}  = naca4((0,0,1,2))`: Foil coordinates at tip
- `position   :: Vector{Real} = zeros(3)`: Position
- `angle      :: Real         = 0.`: Angle of rotation (degrees)
- `axis       :: Vector{Real} = [0.,1.,0.]`: Axis of rotation

### Example

```@example aeromdao
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = Foil.(airfoil for i in 1:3),
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      sweeps      = [0., 30.])

wing = Wing(wing_right);

print_info(wing, "My Wing")
```

You can access each side of a `Wing` by calling either `wing.left` or `wing.right`, and the previous functions should work identically on these `HalfWing`s.

## Fuselage Parametrization