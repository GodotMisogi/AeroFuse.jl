```@meta
EditURL = "<unknown>/docs/lit/howto.jl"
```

# How-to Guide

````@example howto
using AeroMDAO # hide
````

## Airfoil Geometry

### Import Coordinates File
Provide the path to the following function.

````@example howto
# Airfoil coordinates file path
foilpath = string(@__DIR__, "\\..\\..\\data\\airfoil_database\\s1223.dat")

# Read coordinates file
my_foil = read_foil(foilpath)
````

### Interpolate and Process Coordinates

A basic cosine interpolation functionality is provided for foils.

````@example howto
# Cosine spacing with approx. 51 points on upper and lower surfaces each
cos_foil = cosine_spacing(my_foil, 51)
````

The upper and lower surfaces can be obtained by the following variety of functions.

````@example howto
# Upper surface
upper = upper_surface(my_foil)

# Lower surface
lower = lower_surface(my_foil)

# Upper and lower surfaces
upper, lower = split_surface(my_foil)
````

The camber-thickness distribution can be obtained as follows. It internally performs a cosine interpolation to standardize the $x$--coordinates for both surfaces; hence the number of points for the interpolation can be specified.

````@example howto
xcamthick = camber_thickness(my_foil, 60)
````

You can also do the inverse transformation.

````@example howto
coords = camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])
````

## Doublet-Source Aerodynamic Analyses
The `solve_case` method runs the analysis given a `Foil` containing the airfoil coordinates, a `Uniform2D` defining the boundary conditions, and an optional named specification for the number of panels. It returns a system which can be used to obtain the aerodynamic quantities of interest and post-processing.

````@example howto
# Define freestream boundary conditions
uniform = Uniform2D(1.0, 4.0)

# Solve system
system  = @time solve_case(
                     my_foil, uniform;
                     num_panels = 80
                    );

# The following functions compute the quantities of interest, such as the inviscid edge velocities, lift coefficient, and the sectional lift, moment, and pressure coefficients.
panels     = system.surface_panels
@time u_es = surface_velocities(system);
@time cl   = lift_coefficient(system)
@time cls, cms, cps = surface_coefficients(system)
````

AeroMDAO provides more helper functions for the panel geometry.

````@example howto
pts      = collocation_point.(panels) # Collocation point
tangents = panel_tangent.(panels)     # Tangents
normals  = panel_normal.(panels)      # Normals
locs     = panel_location.(panels);   # Upper or lower surface
nothing #hide
````

## Wing Geometry

To define one side of a wing, AeroMDAO provides a `HalfWing` constructor.

````@example howto
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = [ airfoil for i in 1:3 ],
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      LE_sweeps = [0., 30.])
````

The `Wing` constructor takes left and right `HalfWing`s to define a full wing. For example, the following generates a symmetric wing.

````@example howto
wing = Wing(wing_right, wing_right)
````

The following "getter" functions provide quantities of interest such as chord lengths, spans, twist, dihedral, and sweep angles.

```@repl howto
foils(wing)
chords(wing)
spans(wing)
twists(wing)
dihedrals(wing)
sweeps(wing)
```

There is also a convenient function for pretty-printing information, in which the first argument takes the `Wing` type and the second takes a name (as a `String` or `Symbol`).

````@example howto
print_info(wing, "Wing")
````

!!! tip
    You can use [Setfield.jl](https://github.com/jw3126/Setfield.jl) to conveniently copy and modify properties of an existing object.

````@example howto
# Import Setfield
using Setfield

# Set only chords with other properties identical.
wing_left = @set wing_right.chords = [0.4, 0.1, 0.05]
````

To create an asymmetric wing, feed the left and right halves to `Wing` in the particular order.

````@example howto
wing = Wing(wing_left, wing_right);

print_info(wing, "My Wing")
````

## Vortex Lattice Aerodynamic Analyses

How to run a generic aerodynamic analysis on a conventional aircraft configuration.

### Geometry
First we define the lifting surfaces. These can be a combination of `Wing` or `HalfWing` types constructed using the various methods available.
!!! warning "Alert"
    Support for fuselages will be added soon.

````@example howto
# Wing
wing  = WingSection(span       = 8.0,
                    dihedral   = 5.0,
                    LE_sweep   = 15.0,
                    taper      = 0.4,
                    root_chord = 2.0,
                    root_twist = 0.0,
                    tip_twist  = -2.0,
                    root_foil  = naca4(2,4,1,2),
                    tip_foil   = naca4(2,4,1,2),
                    position   = [0., 0., 0.])

# Horizontal tail
htail = WingSection(span       = 2.0,
                    dihedral   = 0.0,
                    LE_sweep   = 15.0,
                    taper      = 0.6,
                    root_chord = 0.8,
                    root_twist = 0.0,
                    tip_twist  = 0.0,
                    root_foil  = naca4(0,0,1,2),
                    tip_foil   = naca4(0,0,0,9),
                    position   = [5., 0., -0.1],
                    angle      = 0.,
                    axis       = [0., 1., 0.]);

# Vertical tail
vtail = HalfWingSection(span       = 0.8,
                        dihedral   = 0.0,
                        LE_sweep   = 8.0,
                        taper      = 0.6,
                        root_chord = 0.8,
                        root_twist = 0.0,
                        tip_twist  = 0.,
                        root_foil  = naca4(0,0,0,9),
                        tip_foil   = naca4(0,0,0,9),
                        position   = [5., 0., 0.],
                        angle      = 90.,
                        axis       = [1., 0., 0.])
````

### Meshing

The `WingMesh` type takes a `HalfWing` or `Wing` type with a vector of integers consisting of the spanwise panel distribution corresponding to the number of sections, and an integer for the chordwise distribution.

````@example howto
# Wing meshes
wing_mesh  = WingMesh(wing, [20], 10)
````

Optionally the type of spanwise spacing can be specified by the keyword `span_spacing` and providing types `Sine(), Cosine(), Uniform()` or a vector of the combination with length corresponding to the number of sections.

````@example howto
# Horizontal tail
htail_mesh = WingMesh(htail, [10], 5;
                      span_spacing = Cosine()) # Uniform(), Sine()

# Vertical tail
vtail_mesh = WingMesh(vtail, [10], 5) # Vertical tail
````

### Inviscid Analysis

The inviscid 3D analysis uses a vortex lattice method. The `WingMesh` type allows you to generate horseshoe elements easily.

````@example howto
wing_horsies = make_horseshoes(wing_mesh)
````

For multiple lifting surfaces, it is most convenient to define a single vector consisting of all the components' horseshoes using [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl).

````@example howto
aircraft = ComponentArray(
                          wing  = wing_horsies,
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );
nothing #hide
````

To define boundary conditions, use the following `Freestream` type, which takes named arguments for angles of attack and sideslip (in degrees), and a quasi-steady rotation vector.

````@example howto
# Define freestream conditions
fs = Freestream(
            alpha = 1.0,
            beta  = 0.0,
            omega = [0.,0.,0.]
           )
````

To define reference values, use the following `References` type.

````@example howto
# Define reference values
refs = References(
           speed     = 10.0,
           density   = 1.225,
           viscosity = 1.5e-5,
           area      = projected_area(wing),
           span      = span(wing),
           chord     = mean_aerodynamic_chord(wing),
           location  = mean_aerodynamic_center(wing)
          )
````

The vortex lattice analysis can be executed with the horseshoes, freestream condition, and reference values defined.

````@example howto
# Solve system
system = solve_case(
             aircraft, fs, refs;
             print            = true, # Prints the results for only the aircraft
             print_components = true, # Prints the results for all components
            )
````

### Dynamics

The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces. You can specify the axis system for the nearfield forces, with choices of geometry, body, wind, or stability. The wind axes are used by default.

````@example howto
# Compute dynamics
ax       = Wind() # Body(), Stability(), Geometry()

# Compute non-dimensional coefficients
CFs, CMs = surface_coefficients(system; axes = ax)
````

You can access the corresponding values of the components' by the name provided in the `ComponentVector`.

````@example howto
CFs.wing
````

Special functions are provided for directly retrieving the dimensionalized forces and moments.

````@example howto
Fs, Ms   = surface_dynamics(system; axes = ax)
````

A Trefftz plane integration is performed to compute farfield forces.

!!! note
    The farfield forces are usually more accurate compared to nearfield forces, as the components do not interact as in the evaluation of the Biot-Savart integrals for the latter.

To obtain the nearfield and farfield coefficients of the components (in wind axes by definition):

````@example howto
nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)
````

You can similarly access the components by name.

````@example howto
@show (nfs.wing, ffs.wing)
````

To obtain the total nearfield and farfield force coefficients:

````@example howto
nf = nearfield(system)
ff = farfield(system)
````

You can also print the relevant information as a pretty table, if necessary.

````@example howto
print_coefficients(nfs.wing, ffs.wing, :wing)
print_coefficients(nf, ff, :aircraft)
````

## Aerodynamic Stability Analyses

````@example howto
ac_name = :aircraft
@time dv_data = solve_case_derivatives(
                       aircraft, fs, refs;
                       axes             = Wind(),
                       name             = ac_name,
                       print            = true,    # Prints the results for only the aircraft
                       print_components = true,    # Prints the results for all components
                      )
````

## Euler-Bernoulli Beam Structural Analysis
The tubular beam's relevant properties, viz. the Young's (elastic) modulus $E$, shear modulus $G$, and torsional moment of inertia $J$ must be specified to define the stiffness matrix for its discretization with $n$ elements.

````@example howto
# Material properties
E = 1. # Elastic modulus
G = 1. # Shear modulus
J = 2. # Torsional moment of inertia
n = 2  # Number of sections

# Stiffness matrix
K = bending_stiffness_matrix(
                             fill(E, 2),
                             fill(G, 2),
                             fill(J, 2),
                             :z         # Direction of deflection
                            )
````

Fixed, hinged beam subjected to force and moment at the center.

````@example howto
# Stiffness matrix
A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3

# Load vector
b = [-1000, 1000, 0]    # F2, M2, M3

# Solution
x = A \ b

# Forces
F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]
````

Propped cantilever beam subjected to force at free end.

````@example howto
# Stiffness matrix
A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2

# Load vector
b = [10, 0, 0]

# Solution
x = A \ b

# Forces
F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

