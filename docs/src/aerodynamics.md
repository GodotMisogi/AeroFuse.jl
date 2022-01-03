# Aerodynamics

```@setup aeromdao
using AeroMDAO
using StaticArrays
using Plots
using LaTeXStrings
gr()
```

The following sections describe different use-cases of AeroMDAO.

## 2D Analyses 

### Airfoil Parametrization

AeroMDAO provides some basic parametrizations commonly used for airfoils, with conversion to the camber-thickness representation and vice versa.

*NACA 4-digit Airfoil Parametrization*

```julia
naca4(digits :: NTuple{4, Real};  # Digits, e.g. (2,4,1,2)
      sharp_trailing_edge = true) # Sharp or blunt trailing edge
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
foil_camthick(coords      :: Array{2, Real}, # 2D coordinates
              num = 40    :: Integer)        # Number of points for distributions 
```

```julia
foilpath  = "path/to/your/airfoil.dat"    # Airfoil coordinates file path
coords    = read_foil(foilpath)           # Read coordinates file
cos_foil  = cosine_foil(coords, 51)       # Cosine spacing with 51 points on upper and lower surfaces
xcamthick = foil_camthick(cos_foil)       # Convert to camber-thickness representation
foiler    = camthick_foil(xcamthick[:,1], # x-components
                          xcamthick[:,2], # Camber distribution
                          xcamthick[:,3]) # Thickness distribution
```

```@example aeromdao
digits  = (2,4,1,2)
airfoil = Foil(naca4(digits))
xcamthick = camber_thickness(airfoil, 60)
coords = camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])
upper, lower = split_foil(coordinates(airfoil));
x_upper, y_upper = upper[:,1], upper[:,2]
x_lower, y_lower = lower[:,1], lower[:,2]
plot(aspectratio = 1)
af_plot = plot(aspect_ratio = 1, xlabel="(x/c)", ylabel = "y")
plot!(x_upper, y_upper, ls = :solid, lw = 2, c = :cornflowerblue, label = "NACA $(digits...) Upper")
plot!(x_lower, y_lower, ls = :solid, lw = 2, c = :orange, label = "NACA $(digits...) Lower")
plot!(xcamthick[:,1], xcamthick[:,2], ls = :dash, lw = 2, c = :burlywood3, label = "NACA $(digits...) Camber")
plot!(xcamthick[:,1], xcamthick[:,3], ls = :dash, lw = 2, c = :grey, label = "NACA $(digits...) Thickness")
```

### Doublet-Source Panel Method

AeroMDAO provides convenience functions using its specific types for analyses.

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

```@example aeromdao
digits  = (0,0,1,2);
airfoil = Foil(naca4(digits; sharp_trailing_edge = true))  #hide
V, α = 1.0, 5.0             #hide
uniform = Uniform2D(V, α)   #hide
cl, cls, cms, cps, panels = solve_case(airfoil, 
                                       uniform;
                                       viscous = false,
                                       sources = false, 
                                       wake_length = 1e3,
                                       wake_panels = 100,
                                       num_panels = 80)

upper, lower = split_foil(coordinates(airfoil)); #hide
x_upper, y_upper = upper[:,1], upper[:,2] #hide
x_lower, y_lower = lower[:,1], lower[:,2] #hide

println("Cl: $cl")
println("Σᵢ Clᵢ: $(sum(cls))")
println("Σᵢ Cmᵢ: $(sum(cms))")
```

```@example aeromdao
digits  = (0,0,1,2); #hide
airfoil = Foil(naca4(digits; sharp_trailing_edge = true))  #hide
V, α = 1.0, 5.0             #hide
uniform = Uniform2D(V, α)   #hide
cl, cls, cms, cps, panels = solve_case(airfoil, #hide
                                       uniform; #hide
                                       viscous = false, #hide
                                       sources = false,  #hide
                                       wake_length = 1e3, #hide
                                       wake_panels = 100, #hide
                                       num_panels = 80) #hide

upper, lower = split_foil(coordinates(airfoil)); #hide
x_upper, y_upper = upper[:,1], upper[:,2] #hide
x_lower, y_lower = lower[:,1], lower[:,2] #hide

get_surface_values(panels, vals, surf = "lower") = [ (collocation_point(panel)[1], val) for (val, panel) in zip(vals, panels) if panel_location(panel) == surf ] #hide
cp_lower = get_surface_values(panels, cps, "lower") #hide
cp_upper = get_surface_values(panels, cps, "upper") #hide
plot(aspect_ratio = 1, yflip = true, xlabel=L"(x/c)", ylabel = L"$C_p$") #hide
plot!(cp_upper, ls = :dash, lw = 2, c = :cornflowerblue, label = "Upper") #hide
plot!(cp_lower, ls = :dash, lw = 2, c = :orange, label = "Lower") #hide
plot!(x_upper, -y_upper,  #hide
        ls = :solid, lw = 2, c = :cornflowerblue, label = "NACA $(digits...) Upper") #hide
plot!(x_lower, -y_lower,  #hide
        ls = :solid, lw = 2, c = :orange, label = "NACA $(digits...) Lower") #hide
```

## 3D Analyses

### Wing Parametrization

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

We can create a `Wing` by feeding two `HalfWing`s to it.

```julia
Wing(left   :: HalfWing, # Left side
     right  :: HalfWing) # Right side
```

```julia
wing_right = HalfWing(foils     = [ Foil(naca4(2,4,1,2)) for i in 1:3 ],
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      LE_sweeps = [0., 30.])

wing       = Wing(wing_right, wing_right)
```

We can obtain the relevant geometric information of the wing for design analyses by calling convenience methods, which automatically perform the necessary calculations on the nested trapezoidal planforms.

```julia
b             = span(wing)
S             = projected_area(wing)
c             = mean_aerodynamic_chord(wing)
AR            = aspect_ratio(wing)
x_w, y_w, z_w = mean_aerodynamic_center(wing)
```

There is also a convenient function for printing this information, whose last optional argument provides the name of the wing. [Pretty Tables](https://github.com/ronisbr/PrettyTables.jl) is used for pretty-printing.

```julia
print_info(wing, "My Wing")
```

```@example aeromdao
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = Foil.(airfoil for i in 1:3),
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      LE_sweeps = [0., 30.])

wing = Wing(wing_right, wing_right);

print_info(wing, "My Wing")
```

You can access each side of a `Wing` by calling either `wing.left` or `wing.right`, and the previous functions should work identically on these `HalfWing`s.


### Vortex Lattice Method

The vortex lattice method used in AeroMDAO follows Mark Drela's *Flight Vehicle Aerodynamics*. The geometry "engine" generates panels for horseshoes and the camber distribution using the airfoil data in the definition of the wing. This geometry is analysed at given freestream angles of attack and sideslip. The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces and a Trefftz plane integration for farfield forces, of which the latter is usually more accurate.

### Meshing
```julia
WingMesh(wing      :: AbstractWing,     # Wing type
         span_num  :: Vector{Integer},  # Number of spanwise panels
         chord_num :: Integer,          # Number of chordwise panels
         spacing   = Cosine())
```

```julia
wing_mesh = WingMesh(wing, [12, 3], 6);
wing_mesh.cam_mesh
```

#### Aircraft Definition

Aircraft analysis by definition of multiple lifting surfaces using the `HalfWing` or `Wing` types are also supported, but is slightly more complex due to the increases in user specifications. Particularly, you will have to mesh the different components by yourself by specifying the number of chordwise panels, spanwise panels and their associated spacings.

This method returns the horseshoe panels for the analysis, and the associated normal vectors based on the camber distribution of the wing. Consider a case in which you have a wing, horizontal tail, and vertical tail.

```julia
# Wing
wing = Wing(foils     = Foil.(fill(naca4(2,4,1,2), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            LE_sweeps = [5.]);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4(0,0,1,2), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4(0,0,0,9), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## Meshing
wing_mesh  = WingMesh(wing, [12], 6, 
                      span_spacing = Cosine()
                     )
htail_mesh = WingMesh(htail, [12], 6, 
                      span_spacing = Cosine()
                     )
vtail_mesh = WingMesh(vtail, [12], 6, 
                      span_spacing = Cosine()
                     )
```

To prepare the analyses, you will have to assemble this information using the exported `ComponentArray`, in which the components are the horseshoes generated using the `make_horseshoes` function.

```julia
aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );
```

Very similarly to the wing-only case, there is an associated `solve_case()` method which takes a dictionary of the panels and normal vectors, and the reference data for computing the non-dimensional coefficients:

#### Freestream Condition

```julia
fs  = Freestream(alpha = 0.0, 
                 beta  = 0.0, 
                 omega = [0., 0., 0.]);
```

##### References

For an analysis, you require ambient reference conditions. Specifically, you need the density $\rho$ and a reference location $r_\text{ref}$ for calculating forces and moments. These dynamics are non-dimensionalized using the wing's characteristics as reference.

```julia
refs = References(speed    = 1.0, 
                  density  = 1.225,
                  area     = projected_area(wing),
                  span     = span(wing),
                  chord    = mean_aerodynamic_chord(wing),
                  location = mean_aerodynamic_center(wing))
```

##### Analysis

```julia
system = solve_case(
                    aircraft, fs, refs;      # Aircraft, Freestream, and References
                    print            = true, # Prints the results for only the aircraft
                    print_components = true, # Prints the results for all components
                   );
```

It returns a system which can be used to determine quantities of interest such as the dynamics."

```julia
# Compute dynamics
ax       = Geometry() # Body(), Stability(), Wind()
CFs, CMs = surface_coefficients(system; axes = ax)
Fs       = surface_forces(system)
Fs, Ms   = surface_dynamics(system; axes = ax) 

nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)

nf  = nearfield(system) 
ff  = farfield(system)
```

##### Stability Analysis

AeroMDAO uses the [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl) package, which leverages forward-mode automatic differentiation to obtain stability derivatives of certain quantities of interest. Particularly, with respect to the angles of attack and sideslip, and the non-dimensionalized roll rates at low computational cost for stability analyses. A unique feature which is not easily obtained from other implementations is the stability derivatives of each component.  To obtain these quantities, simply replace `solve_case()` with `solve_stability_case()`, which will return the data consisting of the nearfield, farfield and stability derivative coefficients.

*Note:* This is due to limitations of using closures in `ForwardDiff`, and hence no post-processing values are provided using this function.

```julia
@time dv_data = solve_stability_case(
                     aircraft, fs, refs;
                     print            = true,    # Prints the results for only the aircraft
                     print_components = true,    # Prints the results for all components
                    );
```

You can pretty-print the stability derivatives with the following function, whose first argument again provides the name of the wing:

```julia
print_derivatives(dv_data.wing, "My Wing")
```



Documentation to be completed. For now, refer to these [analysis](vortex_lattice_method/vlm_aircraft.jl) and [stability analysis](vortex_lattice_method/stability_aircraft.jl) scripts with a full aircraft configuration. There's also an interesting [surrogate model test script](vortex_lattice_method/surrogates.jl)!