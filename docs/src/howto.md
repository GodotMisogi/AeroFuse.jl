```@meta
EditURL = "<unknown>/docs/src/howto.jl"
```

# How-to Guide

````@example howto
using AeroMDAO
````

## Airfoil Processing

### Import Coordinates
Provide the path to the following function.

````@example howto
# Airfoil coordinates file path
foilpath  = string(pwd(), "\\..\\..\\data\\airfoil_database\\s1223.dat")

# Read coordinates file
my_foil    = read_foil(foilpath)
````

### Interpolate Coordinates

````@example howto
# Cosine spacing with 51 points on upper and lower surfaces
cos_foil  = cosine_spacing(my_foil, 51)
````

## Wing Geometry

AeroMDAO provides a `HalfWing` constructor.

````@example howto
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = Foil.(airfoil for i in 1:3),
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      LE_sweeps = [0., 30.])

wing = Wing(wing_right);

print_info(wing, "My Wing")
````

## Vortex Lattice Aerodynamic Analyses

### Meshing

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
# Define meshes
wing_mesh  = WingMesh(wing, [20], 10)
htail_mesh = WingMesh(htail, [10], 5)
vtail_mesh = WingMesh(vtail, [10], 5)
````

### Inviscid Analysis

````@example howto
aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );

# Evaluate case
fs      = Freestream(alpha = 1.0,
                     beta  = 0.0,
                     omega = [0.,0.,0.])
refs    = References(
                     speed    = 10.0,
                     density  = 1.225,
                     area     = projected_area(wing),
                     span     = span(wing),
                     chord    = mean_aerodynamic_chord(wing),
                     location = mean_aerodynamic_center(wing)
                    )

system = solve_case(aircraft, fs, refs;
                    print            = true, # Prints the results for only the aircraft
                    print_components = true, # Prints the results for all components
                   )
````

### Dynamics

The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces. You can specify the axis system for the nearfield forces, with choices of `Geometry(), Body(), Wind(), Stability()`. The wind axes are used by default.

````@example howto
# Compute dynamics
ax       = Wind() # Body(), Stability(), Geometry()
CFs, CMs = surface_coefficients(system; axes = ax)
Fs       = surface_forces(system)
Fs, Ms   = surface_dynamics(system; axes = ax)
````

A Trefftz plane integration is employed to obtain farfield forces.

!!! note
    The farfield forces are usually more accurate compared to nearfield forces, as the components do not interact as in the evaluation of the Biot-Savart integrals for the latter.

To obtain the nearfield and farfield coefficients of the components (in wind axes by definition):

````@example howto
nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)
````

To obtain the total nearfield and farfield force coefficients:

````@example howto
nf  = nearfield(system)
ff  = farfield(system)
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

````@example howto
# Deflection stiffness matrix
K = bending_stiffness_matrix([1., 1.], [1., 1.], [2., 2.], :z)
````

Fixed hinged beam subjected to force and moment at the center.

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

Propped cantilever beam with force at one end.

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

