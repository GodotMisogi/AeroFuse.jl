# Cases

The following sections describe different use-cases of AeroMDAO.

## Airfoil Parametrization

AeroMDAO provides some basic parametrizations commonly used for airfoils, with conversion to the camber-thickness representation and vice versa.

*NACA 4-digit Airfoil Parametrization*

```julia
naca4(digits :: NTuple{4, Real};    # Digits, e.g. (2,4,1,2)
      sharp_trailing_edge = true)   # Sharp or blunt trailing edge
```

```julia
airfoil = naca4((2,4,1,2))
```

*Kulfan Class Shape Transformation (CST) Method*
```julia
kulfan_CST(alpha_u      :: Vector{Real},    # Upper surface parameters
           alpha_l      :: Vector{Real},    # Lower surface parameters
           dzs          :: NTuple{2, Real}, # Upper and lower trailing edge points
           coeff_LE = 0 :: Real,            # Leading-edge modification coefficient
           n = 40       :: Integer)         # Number of points on each surface
```

```julia
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
dzs     = (1e-4, 1e-4)
foil    = kulfan_CST(alpha_u, alpha_l, dzs, 0.2)
```

*Camber-Thickness Representation*
```julia
foil_camthick(coords	  :: Array{2, Real},  # 2D coordinates
              num = 40	  :: Integer)       # Number of points for distributions
```

```julia
foilpath  = "path/to/your/airfoil.dat"       # Airfoil coordinates file path
coords    = read_foil(foilpath)              # Read coordinates file
cos_foil  = cosine_foil(coords, 51)          # Cosine spacing with 51 points on upper and lower surfaces
xcamthick = camber_thickness(Foil(cos_foil)) # Convert to camber-thickness representation
foiler    = camber_thickness_to_coordinates(xcamthick[:,1], # x-components
                          xcamthick[:,2],    # Camber distribution
                          xcamthick[:,3])    # Thickness distribution
```

## Doublet-Source Panel Method

AeroMDAO provides convenient functions using its specific types for analyses.

The `Foil` type converts the coordinates into a friendly type for analyses.

```julia
airfoil = naca4((2,4,1,2))
foil    = Foil(airfoil)
```

The `Uniform2D` type consists of the freestream speed and angle of attack.

```julia
V, α    = 1.0, 3.0
uniform = Uniform2D(V, α)
```

To analyse this foil with these boundary conditions using the incompressible 2D doublet-source panel method, the following method is called. Optional named arguments are provided to specify whether the source terms are non-zero, the length of the wake, and the number of panels for the analysis.

```julia
solve_case(foil            :: Foil,
           uniform         :: Uniform2D;
           num_panels = 60 :: Integer)
```

The method returns the lift coefficient calculated by the doublet strength of the wake panel, the lift, moment and pressure coefficients over the panels, and the panels generated for post-processing.

```julia
cl, cls, cms, cps, panels =
  solve_case(foil,
             uniform;
             num_panels = 80)
```

_Disclaimer_: Only sharp trailing edge airfoils are supported for now. Blunt trailing edge support will be added soon.

## Wing Parametrization

The following image depicts the parametrization schema used for wing planforms in terms of nested trapezoids.

![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

This parametrization is implemented as a composite type called a `HalfWing`.

```julia
HalfWing(foils     :: Vector{Foil}, # Foil profiles
         chords	   :: Vector{Real}, # Chord lengths
         twists	   :: Vector{Real}, # Twist angles (deg)
         spans	   :: Vector{Real}, # Section span lengths
         dihedrals :: Vector{Real}, # Dihedral angles (deg)
         sweeps	   :: Vector{Real}, # Leading-edge sweep angles (deg)
         position  = [0., 0., 0.],  # Global leading edge position
         angle     = 0.,            # Angle of rotation (deg)
         axis      = [0., 1., 0.],  # Axis/plane of rotation about leading edge
```

We can create a `Wing` by feeding two `HalfWing`s to it.

```julia
Wing(left   :: HalfWing, # Left side
     right  :: HalfWing) # Right side
```

```julia
airfoil    = naca4((2,4,1,2))
foils      = Foil.(airfoil for i in 1:3)
wing_right = HalfWing(foils,
                      [0.4, 0.2, 0.1],
                      [0., 2., 5.],
                      [1.0, 0.1],
                      [0., 60.],
                      [0., 30.])

wing       = Wing(wing_right, wing_right)
```

We can obtain the relevant geometric information of the wing for design analyses by calling convenience methods, which automatically perform the necessary calculations on the nested trapezoidal planforms.

```julia
b   = span(wing)
S   = projected_area(wing)
AR  = aspect_ratio(wing)
c   = mean_aerodynamic_chord(wing)
mac = mean_aerodynamic_center(wing)
```

There is also a convenient function for printing this information, whose last optional argument provides the name of the wing. [Pretty Tables](https://github.com/ronisbr/PrettyTables.jl) is used for pretty-printing.

```julia
print_info(wing, "My Wing")
```

You can access each side of a `Wing` by calling either `wing.left` or `wing.right`, and the previous functions should work identically on these `HalfWing`s.


## Vortex Lattice Method

The vortex lattice method used in AeroMDAO follows Mark Drela's *Flight Vehicle Aerodynamics*. The geometry "engine" generates panels for horseshoes and the camber distribution using the airfoil data in the definition of the wing. This geometry is analysed at given freestream angles of attack and sideslip. The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces and a Trefftz plane integration for farfield forces, of which the latter is usually more accurate.

### Wing Analysis

To run a vortex lattice analysis on a `Wing`, we first define a `Freestream` object consisting of the boundary conditions for the analyses, parametrised by the freestream speed `V`, angle of attack `α`, side-slip angle `β`, and a quasi-steady rotation vector `Ω`. We also specify the density `ρ` and reference location `r_ref` for calculating moments generated by the forces.

```julia
U = 10.0
α = 5.0
β = 0.0
Ω = [0.0, 0.0, 0.0]

freestream = Freestream(U, α, β, Ω);

ρ     = 1.225
r_ref = [ mean_aerodynamic_center(wing)[1], 0., 0.]
```

Now we run the case with specifications of the number of spanwise and chordwise panels by calling the `solve_case()` function, which has an associated method.

```julia
solve_case(wing                   :: Union{Wing, HalfWing},
           freestream             :: Freestream;
           r_ref = [0.25, 0., 0.] :: Vector{Real},
           rho_ref = 1.225        :: Real,
           area_ref  = 1.,        :: Real,
           span_ref  = 1.,        :: Real,
           chord_ref = 1.,        :: Real,
           span_num = 25          :: Union{Integer, Vector{Integer}},
           chord_num = 10         :: Integer,
           viscous = false        :: Bool,
           x_tr = 0.3             :: Union{Real, Vector{Real}}
          )
```

The method uses the planform of the wing to compute the horseshoe elements, and the camber distribution to compute the normal vectors for the analyses. The user can specify the distribution of spanwise panels as a vector depending on the number of sections, and the relevant reference values. If geometric ones are not provided, i.e. the reference span, chord, and area, then the function will automatically calculate them using the properties of the wing as shown previously.

A "viscous" analysis is also supported using traditional wetted-area methods based on the equivalent skin-friction drag formulation of Schlichting, which have limited applicability due to their reliance on empirical form factors. The transition locations over each section can be specified as a vector of fractions of the normalized local average chord lengths of each section, or alternatively as a number for fixing it at approximately the same local average fraction over all sections.

The analysis returns the nearfield and farfield coefficients, the non-dimensionalised force and moment coefficients over the panels, the panels used for the horseshoe elements, the camber panels, the horseshoe elements used for plotting streamlines, and the associated vortex strengths `Γs` corresponding to the solution of the system. The following [example script](vortex_lattice_method/vlm_wing.jl) is provided for reference.
```julia
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normies, horses, Γs =
    solve_case(wing, freestream;
               rho_ref   = ρ,
               r_ref     = r_ref,
               area_ref  = S,
               span_ref  = b,
               chord_ref = c,
               span_num  = [15, 9],
               chord_num = 6,
               viscous   = true,
               x_tr      = [0.3, 0.3]);
```

You can pretty-print the aerodynamic coefficients with the following function, whose last optional argument provides the name of the wing:

```julia
print_coefficients(nf_coeffs, ff_coeffs, "My Wing")
```

If the viscous option is enabled, this function also prints the pressure, induced, and total drags separately.

#### Stability Derivatives

AeroMDAO uses the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package, which leverages forward-mode automatic differentiation to obtain the stability derviatives with respect to the angle of attack and sideslip and the non-dimensionalized roll rates at low computational cost. To obtain the stability derivatives, simply replace `solve_case()` with `solve_stability_case()`, which will return **only** the nearfield, farfield and stability derivative coefficients. This is due to limitations of using closures in `ForwardDiff`, and hence no post-processing values are provided using this function.

```julia
nf_coeffs, ff_coeffs, dv_coeffs =
    solve_stability_case(wing, freestream;
                         rho_ref   = ρ,
                         r_ref     = r_ref,
                         area_ref  = S,
                         span_ref  = b,
                         chord_ref = c,
                         span_num  = [15, 9],
                         chord_num = 6,
                         viscous   = false,
                         x_tr      = [0.3, 0.3]);
```

You can pretty-print the stability derivatives with the following function, whose last argument again optioanlly provides the name of the wing:

```julia
print_derivatives(dv_coeffs, "Wing")
```

### Aircraft Analysis

Aircraft analysis by definition of multiple lifting surfaces using the `HalfWing` or `Wing` types are also supported, but is slightly more complex due to the increases in user specifications. Particularly, you will have to mesh the different components yourself by specifying the number of chordwise panels, spanwise panels and their associated spacing distributions.

```julia
panel_wing(wing                :: Union{HalfWing, Wing},
           span_num,           :: Union{Integer, Vector{<: Integer}},
           chord_num           :: Integer;
           spacing  = "cosine" :: Union{String, Vector{String}} # Options: "uniform", "sine", "cosine"
          )
```

This method returns the horseshoe panels for the analysis, and the associated normal vectors based on the camber distribution of the wing. Consider a case in which you have a wing, horizontal tail, and vertical tail.

```julia
# Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 3)),
            chords    = [1.0, 0.6, 0.2],
            twists    = [0.0, 0.0, 0.0],
            spans     = [5.0, 0.5],
            dihedrals = [5., 5.],
            LE_sweeps = [5., 5.]);

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
htail = Wing(foils     = htail_foils,
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail_foils = Foil.(fill(naca4((0,0,0,9)), 2))
vtail = HalfWing(foils     = vtail_foils,
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

# Assembly
wing_panels , wing_normals  = panel_wing(wing, [20, 5], 10;
                                         spacing = "uniform"
                                        )
htail_panels, htail_normals = panel_wing(htail, [6], 6);
                                         spacing = "uniform"
                                        )
vtail_panels, vtail_normals = panel_wing(vtail, [6], 5;
                                         spacing = "uniform"
                                        )
```

To prepare the analyses, you will have to assemble this information into a dictionary, where the keys are the names of the components.

```julia
aircraft = Dict("Wing"            => Horseshoe.(wing_panels,  wing_normals ),
                "Horizontal Tail" => Horseshoe.(htail_panels, htail_normals),
                "Vertical Tail"   => Horseshoe.(vtail_panels, vtail_normals))
```

Very similarly to the wing-only case, there is an associated `solve_case()` method which takes a dictionary of the panels and normal vectors, and the reference data for computing the non-dimensional coefficients:

```julia
ac_name = "My Aircraft"
ρ       = 1.225
V, α, β = 1.0, 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)

S, b, c  = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
wing_mac = mean_aerodynamic_chord(wing) 
x_w      = wing_mac[1]
ref      = [ x_w, 0., 0.]

data =
    solve_case(aircraft, fs;
               rho_ref          = ρ,        # Reference density
               r_ref            = ref,      # Reference point for moments
               area_ref         = S,        # Reference area
               span_ref         = b,        # Reference span
               chord_ref        = c,        # Reference chord
               name             = ac_name,  # Aircraft name
               print            = true,     # Prints the results for the entire aircraft
               print_components = true,     # Prints the results for each component
              );
```

And similarly for obtaining the stability derivatives. A unique feature which is not easily obtained from other implementations (really? should check...) is the stability derivatives of each component thanks to the automatic differentiation setup:

```julia
dv_data =
    solve_stability_case(aircraft, fs;
                         rho_ref          = ρ,
                         r_ref            = ref,
                         area_ref         = S,
                         span_ref         = b,
                         chord_ref        = c,
                         name             = ac_name,
                         print            = true,
                         print_components = true,
                        );
```

Documentation to be completed. For now, refer to these [analysis](vortex_lattice_method/vlm_aircraft.jl) and [stability analysis](vortex_lattice_method/stability_aircraft.jl) scripts with a full aircraft configuration. There's also an interesting test [surrogate model test script](vortex_lattice_method/surrogates.jl)!

## Aerostructural Analysis
