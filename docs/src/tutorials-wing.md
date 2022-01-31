```@meta
EditURL = "<unknown>/docs/lit/tutorials-wing.jl"
```

## Objectives

Here we will teach you the basic functionality of AeroMDAO by showing you how to perform an aerodynamic analysis of a conventional aircraft. For this, we will need to import some packages which will be convenient for plotting.

````@example tutorials-wing
using AeroMDAO      # Main package
using Plots         # Plotting library
gr(dpi = 300)       # Plotting backend
using LaTeXStrings  # For LaTeX printing in plots
````

## Your First Wing

Here you will learn how to define a wing using an intuitive parametrization scheme. First, we define a `Vector` of `Foil`s.

````@example tutorials-wing
# Define one airfoil
airfoil_1 = naca4((4,4,1,2))

# Define vector of airfoils
airfoils  = [ airfoil_1, naca4((0,0,1,2)) ]
````

!!! note
    Refer to the [Airfoil Aerodynamic Analysis](tutorials-airfoil.md) tutorial for an introduction to the `Foil` type.

### Parametrization

The following parametrization is used for the wing, presented for visual understanding.
![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

The following function defines a symmetric wing and prints the relevant information.

````@example tutorials-wing
wing = Wing(foils     = airfoils,    # Foil profiles
            chords    = [1.0, 0.6],  # Chord lengths
            twists    = [2.0, 0.0],  # Twist angles (degrees)
            spans     = [4.0],       # Section span lengths
            dihedrals = [5.],        # Dihedral angles (degrees)
            LE_sweeps = [5.])        # Leading-edge sweep angles (degrees)
````

!!! note
    See the [how-to guide](howto.md) on how to define an asymmetric wing.

### Visualization

The following function generates the coordinates of the wing's outline.

````@example tutorials-wing
wing_outline = plot_wing(wing)
````

Let's plot the geometry!

````@example tutorials-wing
plt = plot(
           wing_outline[:,1], wing_outline[:,2], wing_outline[:,3],
           label = "Wing",
           xaxis = "x", yaxis = "y", zaxis = "z",
           aspect_ratio = 1,
           camera = (30, 45),
           zlim = (-0.1, span(wing) / 2),
          )
````

## Your First Vortex Lattice Analysis
Now we would like to analyze the aerodynamics of this wing in conjunction with other lifting surfaces pertinent to aircraft.

### Geometry
We define the horizontal tail similarly to the wing.

````@example tutorials-wing
# Horizontal tail
htail = Wing(foils     = fill(naca4(0,0,1,2), 2),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             LE_sweeps = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])
````

For the vertical tail, we simply replace `Wing` with `HalfWing` to define its shape.

````@example tutorials-wing
# Vertical tail
vtail = HalfWing(foils     = fill(naca4(0,0,0,9), 2),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])
````

Let's visualize the geometry of the aircraft's configuration.

````@example tutorials-wing
htail_outline = plot_wing(htail)
vtail_outline = plot_wing(vtail)

plot!(plt,
    htail_outline[:,1], htail_outline[:,2], htail_outline[:,3],
    label = "Horizontal Tail"
   )
plot!(plt,
   vtail_outline[:,1], vtail_outline[:,2], vtail_outline[:,3],
   label = "Vertical Tail"
  )
````

### Meshing and Assembly
To perform the aerodynamic analysis, we will need to discretize our geometry into a _mesh_. The following `WingMesh` function constructs a mesh for you by providing a `HalfWing` or `Wing` type with specification of spanwise panels and chordwise panels. As the wing has only one spanwise section, we provide a vector with a single  integer entry for the spanwise panel distribution.

````@example tutorials-wing
wing_mesh = WingMesh(wing, [12], 6) # (Wing, [Spanwise panels], Chordwise panels)
````

Let's see what this discretization looks like on the camber distribution of the wing.

````@example tutorials-wing
# Compute camber panel distribution
wing_cam_panels = camber_panels(wing_mesh)

# Generate plotting points
plt_wing_pans   = plot_panels(wing_cam_panels)

# Plot panels
plt_pans = plot(
             xaxis = "x", yaxis = "y", zaxis = "z",
             aspect_ratio = 1,
             camera = (30, 45),
             zlim = (-0.1, span(wing) / 2),
            )

[ plot!(plt_pans, panel, label = "", color = :lightblue) for panel in plt_wing_pans ]
plot!(plt_pans)
````

Similarly we define the meshes for the other surfaces and plot them.

````@example tutorials-wing
htail_mesh = WingMesh(htail, [12], 6)
vtail_mesh = WingMesh(vtail, [12], 6)

[ plot!(plt_pans, panel, label = "", color = :orange) for panel in plot_panels(camber_panels(htail_mesh)) ]
[ plot!(plt_pans, panel, label = "", color = :lightgreen) for panel in plot_panels(camber_panels(vtail_mesh)) ]
plot!(plt_pans)
````

For the analysis, you have to assemble the meshes.

````@example tutorials-wing
aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         )
````

Define freestream condition.

````@example tutorials-wing
fs  = Freestream(
                 alpha = 3.0, # degrees
                 beta  = 0.0, # degrees
                 omega = [0., 0., 0.]
                );
nothing #hide
````

Define reference values.

````@example tutorials-wing
refs = References(speed    = 1.0,
                  area     = projected_area(wing),
                  span     = span(wing),
                  chord    = mean_aerodynamic_chord(wing),
                  density  = 1.225,
                  location = mean_aerodynamic_center(wing));
nothing #hide
````

The

````@example tutorials-wing
system = solve_case(
        aircraft, fs, refs;
        print            = true, # Prints the results for only the aircraft
        print_components = true, # Prints the results for all components
    )
````

Y

````@example tutorials-wing
dv_data = solve_case_derivatives(
        aircraft, fs, refs;
        print              = false, # Prints the results for only the aircraft
        print_components   = true,  # Prints the results for all components
    );
nothing #hide
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

