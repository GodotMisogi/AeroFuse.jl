# # How-to Guide
using AeroMDAO # hide
using Plots # hide
gr(dpi = 300) # hide
using LaTeXStrings # hide

# ## Airfoil Geometry
# How to work with airfoil geometry.

# ### Import Coordinates File
# You can specify the path consisting of the foil's coordinates to the `read_foil` function. The format requires a header for the name by default, but this can be disabled and a custom name can be provided by setting the optional `header` and `name` variables.

## Airfoil coordinates file path
foilpath = string(@__DIR__, "/misc/s1223.dat")

## Read coordinates file
my_foil = read_foil(foilpath;
                    header = true,
                    name   = "")

#
plot(my_foil.x, my_foil.y, 
     xlabel = L"x", ylabel = L"y", aspect_ratio = 1, label = "$(my_foil.name)")

# ### Interpolate and Process Coordinates

# A basic cosine interpolation functionality is provided for foils.
## Cosine spacing with approx. 51 points on upper and lower surfaces each
cos_foil = cosine_interpolation(my_foil, 51)

# The upper and lower surfaces can be obtained by the following variety of functions.

## Upper surface
upper = upper_surface(my_foil)

## Lower surface
lower = lower_surface(my_foil)

## Upper and lower surfaces
upper, lower = split_surface(my_foil)

# The camber-thickness distribution can be obtained as follows. It internally performs a cosine interpolation to standardize the $x$--coordinates for both surfaces; hence the number of points for the interpolation can be specified.
xcamthick = camber_thickness(my_foil, 60)

# You can also do the inverse transformation.
coords = camber_thickness_to_coordinates(xcamthick[:,1], xcamthick[:,2], xcamthick[:,3])

# ### Control Surfaces

# You can (somewhat) mimic the behaviour of a control surface by specifying a deflection angle $\delta$ (in degrees, clockwise-positive convention) with the specification of the hinge location's $x$-coordinate normalized in $[0,1]$ to the chord length.
con_foil = control_surface(cos_foil; angle = 10., hinge = 0.75)

plot!(con_foil.x, con_foil.y, label = "$(my_foil.name) Deflected")

# ## Doublet-Source Aerodynamic Analyses
# The `solve_case` method runs the analysis given a `Foil` containing the airfoil coordinates, a `Uniform2D` defining the boundary conditions, and an optional named specification for the number of panels. It returns a system which can be used to obtain the aerodynamic quantities of interest and post-processing.

## Define freestream boundary conditions
uniform = Uniform2D(1.0, 4.0)

## Solve system
system  = solve_case(my_foil, uniform;
                     num_panels = 80)

# The following functions compute the quantities of interest, such as the inviscid edge velocities, lift coefficient, and the sectional lift, moment, and pressure coefficients.
cls, cms, cps = surface_coefficients(system);
u_es   = surface_velocities(system)
cl     = lift_coefficient(system)

# AeroMDAO provides more helper functions for the panel geometry.
panels   = system.surface_panels
pts      = collocation_point.(panels) # Collocation points
tangents = panel_tangent.(panels)     # Tangents
normals  = normal_vector.(panels)      # Normals
locs     = panel_location.(panels);   # Upper or lower surface

# ## Wing Geometry
# How to work with wing geometry.
# 
# To define one side of a wing, AeroMDAO provides a `HalfWing` constructor.
airfoil    = naca4((2,4,1,2))
wing_right = HalfWing(foils     = [ airfoil for i in 1:3 ],
                      chords    = [0.4, 0.2, 0.1],
                      twists    = [0., 2., 5.],
                      spans     = [1.0, 0.1],
                      dihedrals = [0., 60.],
                      sweeps    = [0., 30.],
                      w_sweep   = 0.25)

# The `Wing` constructor takes left and right `HalfWing`s to define a full wing. For example, the following generates a symmetric wing.
wing = Wing(wing_right, wing_right)

# The following "getter" functions provide quantities of interest such as chord lengths, spans, twist, dihedral, and sweep angles.

# ```@repl howto
# foils(wing)
# chords(wing)
# spans(wing)
# twists(wing)
# dihedrals(wing)
# sweeps(wing)
# ```

# There is also a convenient function for pretty-printing information, in which the first argument takes the `Wing` type and the second takes a name (as a `String` or `Symbol`).
print_info(wing, "Wing")

#md # !!! tip
#md #     You can use [Setfield.jl](https://github.com/jw3126/Setfield.jl) to conveniently copy and modify properties of an existing object.

## Import Setfield
using Setfield

## Set only chords with other properties remaining identical.
wing_left = @set wing_right.chords = [0.4, 0.1, 0.05]

# To create an asymmetric wing, feed the left and right halves to `Wing` in the particular order.
wing = Wing(wing_left, wing_right);

print_info(wing, "My Wing")

# ## Vortex Lattice Aerodynamic Analyses
# 
# How to run a generic aerodynamic analysis on a conventional aircraft configuration.
#
# ### Geometry 
# First we define the lifting surfaces. These can be a combination of `Wing` or `HalfWing` types constructed using the various methods available.
# !!! warning "Alert"
#     Support for fuselages and control surfaces will be added soon.

## Wing
wing  = WingSection(span       = 8.0,
                    dihedral   = 5.0,
                    sweep      = 15.0,
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
                    sweep      = 15.0,
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
                        sweep      = 8.0,
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
vtail_mesh = WingMesh(vtail, [10], 5); # Vertical tail

# ### Inviscid Analysis

# The inviscid 3D analysis uses a vortex lattice method. The `WingMesh` type allows you to generate horseshoe elements easily.
wing_horsies = make_horseshoes(wing_mesh)

# For multiple lifting surfaces, it is most convenient to define a single vector consisting of all the components' horseshoes using [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl).
aircraft = ComponentVector(
                           wing  = wing_horsies,
                           htail = make_horseshoes(htail_mesh),
                           vtail = make_horseshoes(vtail_mesh)
                          );

# To define boundary conditions, use the following `Freestream` type, which takes named arguments for angles of attack and sideslip (in degrees), and a quasi-steady rotation vector.

## Define freestream conditions
fs = Freestream(alpha = 1.0, 
                beta  = 0.0, 
                omega = [0.,0.,0.])

# To define reference values, use the following `References` type.

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

# The vortex lattice analysis can be executed with the horseshoes, freestream condition, and reference values defined.

## Solve system
system = solve_case(
             aircraft, fs, refs;
             print            = true, # Prints the results for only the aircraft
             print_components = true, # Prints the results for all components
            )

# If needed, you can access the relevant component influence matrix values and boundary conditions with the following attributes, e.g.
system.influence_matrix[:wing] # Or :htail, :vtail
system.influence_matrix[:wing, :htail];

system.boundary_vector[:wing]

# ### Dynamics
# 
# The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces. You can specify the axis system for the nearfield forces, with choices of geometry, body, wind, or stability axes. The wind axes are used by default.

## Compute dynamics
ax = Wind() # Body(), Stability(), Geometry()

## Compute non-dimensional coefficients
CFs, CMs = surface_coefficients(system; axes = ax)

# You can access the corresponding values of the components' by the name provided in the `ComponentVector`.
CFs.wing

# Special functions are provided for directly retrieving the dimensionalized forces and moments.
Fs, Ms   = surface_dynamics(system; axes = ax)

# A Trefftz plane integration is performed to compute farfield forces.
# 
# !!! note
#     The farfield forces are usually more accurate compared to nearfield forces, as the components do not interact as in the evaluation of the Biot-Savart integrals for the latter.
# 
# To obtain the nearfield and farfield coefficients of the components (in wind axes by definition):
nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)

# You can similarly access the components by name.
@show (nfs.wing, ffs.wing)

# To obtain the total nearfield and farfield force coefficients:
nf = nearfield(system) 
ff = farfield(system)

# You can also print the relevant information as a pretty table, if necessary.
print_coefficients(nfs.wing, ffs.wing, :wing)
print_coefficients(nf, ff, :aircraft)

# ## Aerodynamic Stability Analyses
# The derivatives of the aerodynamic coefficients with respect to the freestream values is obtained by automatic differentiation enabled by [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl). To compute the values, simply replace `solve_case` with `solve_case_derivatives`. You can also optionally provide the axes for the reference frame of the coefficients.
dv_data = solve_case_derivatives(
             aircraft, fs, refs;
             axes             = Wind(),
             name             = :aircraft,
             print            = true,    # Prints the results for only the aircraft
             print_components = true,    # Prints the results for all components
            );

# ## Euler-Bernoulli Beam Structural Analysis
# The tubular beam's relevant properties, viz. the Young's (elastic) modulus $E$, shear modulus $G$, and torsional moment of inertia $J$ must be specified to define the stiffness matrix for its discretization with $n$ elements.

## Material properties
E = 1. # Elastic modulus
G = 1. # Shear modulus
J = 2. # Torsional moment of inertia
n = 2  # Number of sections

## Stiffness matrix
K = bending_stiffness_matrix(
                             fill(E, 2), 
                             fill(G, 2),
                             fill(J, 2), 
                             :z         # Direction of deflection
                            )

# Fixed, hinged beam subjected to force and moment at the center.

## Stiffness matrix
A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3

## Load vector
b = [-1000, 1000, 0]    # F2, M2, M3

## Solution
x = A \ b

## Forces
F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]

# Propped cantilever beam subjected to force at free end.

## Stiffness matrix
A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2

## Load vector
b = [10, 0, 0]

## Solution
x = A \ b

## Forces
F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]
