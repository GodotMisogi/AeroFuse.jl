# # How-to Guide
using AeroMDAO # hide

# ## Airfoil Processing

# ### Import Coordinates
# Provide the path to the following function.

## Airfoil coordinates file path
foilpath = string(@__DIR__, "\\..\\..\\data\\airfoil_database\\s1223.dat")

## Read coordinates file
my_foil = read_foil(foilpath)

# ### Interpolate and Process Coordinates

# A basic cosine interpolation functionality is provided for foils.
## Cosine spacing with approx. 51 points on upper and lower surfaces each
cos_foil = cosine_spacing(my_foil, 51)

# The upper and lower surfaces can be obtained by the following variety of functions.

## Upper surface
upper = upper_surface(my_foil)

## Lower surface
lower = lower_surface(my_foil)

## Upper and lower surfaces
upper, lower = split_surface(my_foil)

# The camber-thickness distribution can be obtained as follows. It internally performs a cosine interpolation to standardize the `$x$--coordinates for both surfaces; hence the number of points for the interpolation can be specified.
xcamthick = camber_thickness(my_foil, 60)

# ## Wing Geometry
# 
# AeroMDAO provides a `HalfWing` constructor .
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = Foil.(airfoil for i in 1:3),
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      LE_sweeps = [0., 30.])

#md # !!! tip
#md #     You can use [Setfield.jl](https://github.com/jw3126/Setfield.jl) to conveniently copy and modify properties.

##
using Setfield

## 
wing_left = @set wing_right.chords = [0.4, 0.1, 0.05]

# To create an asymmetric wing, feed the left and right halves to `Wing` in the particular order.
wing = Wing(wing_left, wing_right);

print_info(wing, "My Wing")

# ## Doublet-Source Aerodynamic Analyses
# The method returns the lift coefficient calculated by the doublet strength of the wake panel, the lift, moment and pressure coefficients over the panels, and the panels generated for post-processing.
cl, cls, cms, cps, panels = solve_case(foil, uniform; num_panels = 80)

# ## Vortex Lattice Aerodynamic Analyses
# 
# How to run a generic aerodynamic analysis on a conventional aircraft configuration.
#
# ### Geometry 

## Wing
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

## Horizontal tail
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

## Vertical tail
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

# ### Meshing

# The `WingMesh` type takes a `HalfWing` or `Wing` type with a vector of integers consisting of the spanwise panel distribution corresponding to the number of sections, and an integer for the chordwise distribution.

## Wing meshes
wing_mesh  = WingMesh(wing, [20], 10)

# Optionally the type of spanwise spacing can be specified by the keyword `span_spacing` and providing types `Sine(), Cosine(), Uniform()` or a vector of the combination with length corresponding to the number of sections.

## Horizontal tail
htail_mesh = WingMesh(htail, [10], 5;
                      span_spacing = Cosine()) # Uniform(), Sine()

## Vertical tail
vtail_mesh = WingMesh(vtail, [10], 5) # Vertical tail

# ### Inviscid Analysis

# The inviscid 3D analysis uses a vortex lattice method. The `WingMesh` type allows you to generate horseshoe elements easily.
wing_horsies = make_horseshoes(wing_mesh)

# For multiple lifting surfaces, it is most convenient to define a single vector consisting of all the components' horseshoes using [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl).
aircraft = ComponentArray(
                          wing  = wing_horsies,
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );

# To define boundary conditions, use the following `Freestream` type, which takes named arguments for angles of attack and sideslip (in degrees), and a quasi-steady rotation vector.

## Define freestream conditions
fs = Freestream(
            alpha = 1.0, 
            beta  = 0.0, 
            omega = [0.,0.,0.]
           )

# To define reference values, use the following `References` type 

## Define reference values
refs = References(
           speed     = 10.0,
           density   = 1.225,
           viscosity = 1.5e-5,
           area      = projected_area(wing),
           span      = span(wing), 
           chord     = mean_aerodynamic_chord(wing), 
           location  = mean_aerodynamic_center(wing)
          )

## Solve system
system = solve_case(
             aircraft, fs, refs;
             print            = true, # Prints the results for only the aircraft
             print_components = true, # Prints the results for all components
            )

# ### Dynamics
# 
# The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces. You can specify the axis system for the nearfield forces, with choices of geometry, body, wind, or stability. The wind axes are used by default.
# 

## Compute dynamics
ax       = Wind() # Body(), Stability(), Geometry()

## Compute non-dimensional coefficients
CFs, CMs = surface_coefficients(system; axes = ax)

# Special functions are provided for directly retrieving the dimensionalized forces and moments.
Fs, Ms   = surface_dynamics(system; axes = ax)

# A Trefftz plane integration is employed to obtain farfield forces.
# 
# !!! note
#     The farfield forces are usually more accurate compared to nearfield forces, as the components do not interact as in the evaluation of the Biot-Savart integrals for the latter.
# 
# To obtain the nearfield and farfield coefficients of the components (in wind axes by definition):
nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)

# To obtain the total nearfield and farfield force coefficients:
nf  = nearfield(system) 
ff  = farfield(system)

# ## Aerodynamic Stability Analyses

ac_name = :aircraft
@time dv_data = solve_case_derivatives(
                       aircraft, fs, refs;
                       axes             = Wind(),
                       name             = ac_name,
                       print            = true,    # Prints the results for only the aircraft
                       print_components = true,    # Prints the results for all components
                      )


# ## Euler-Bernoulli Beam Structural Analysis

## Deflection stiffness matrix
K = bending_stiffness_matrix([1., 1.], [1., 1.], [2., 2.], :z)

# Fixed hinged beam subjected to force and moment at the center.

## Stiffness matrix
A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3

## Load vector
b = [-1000, 1000, 0]    # F2, M2, M3

## Solution
x = A \ b

## Forces
F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]

# Propped cantilever beam with force at one end.

## Stiffness matrix
A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2

## Load vector
b = [10, 0, 0]

## Solution
x = A \ b

## Forces
F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]
