# # How-to Guide
using AeroFuse # hide
using Plots # hide
gr(dpi = 300) # hide
using LaTeXStrings # hide

# ## Airfoil Geometry
# How to work with airfoil geometry.

# ### Import Coordinates File
# You can specify the path consisting of the foil's coordinates to the `read_foil` function. The [Selig format](https://openvsp.org/wiki/doku.php?id=airfoilexport) for the coordinates file is followed, in which a header for the name in the first line of the file with the coordinates following from the second line are required. A custom name can be provided by setting the optional `name` variable.

## Airfoil coordinates file path
foilpath = string(@__DIR__, "/misc/s1223.dat")

# Note that, in this example, `@__DIR__` points to the Julia REPL's current working directory. Use the correct local path for the coordinates on your computer.

## Read coordinates file
my_foil = read_foil(
    foilpath; # Path to the airfoil coordinates file
    name = "S1223" # To overwrite name in header
)

#
plot(xlabel = L"x", ylabel = L"y", aspect_ratio = 1)
plot!(my_foil)

# !!! tip
#     You can also import airfoil coordinates from the internet.

## Read coordinates from an internet address.
online_foil = read_foil(download("https://m-selig.ae.illinois.edu/ads/coord/s1210.dat"))

plot!(online_foil)

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

plot!(con_foil)

# ## Doublet-Source Aerodynamic Analyses
# The `solve_case` method runs the analysis given a `Foil` containing the airfoil coordinates, a `Uniform2D` defining the boundary conditions, and an optional named specification for the number of panels. It returns a system which can be used to obtain the aerodynamic quantities of interest and post-processing.

## Define freestream boundary conditions
uniform = Uniform2D(1.0, 4.0)

## Solve system
system  = solve_case(
    my_foil, # Foil
    uniform; # Freestream condition
    num_panels = 80
)

# The following functions compute the quantities of interest, such as the inviscid edge velocities, lift coefficient, and the sectional lift, moment, and pressure coefficients.
cls, cms, cps = surface_coefficients(system);
u_es   = surface_velocities(system)
cl     = lift_coefficient(system)

# AeroFuse provides more helper functions for the panel geometry.
panels   = system.surface_panels
pts      = collocation_point.(panels) # Collocation points
tangents = tangent_vector.(panels)     # Tangents
normals  = normal_vector.(panels)      # Normals
locs     = panel_location.(panels);   # Upper or lower surface

# ## Wing Geometry
# How to work with wing geometry.
# 
# To define a wing, AeroFuse provides a `Wing` constructor based on the following parametrization. The named arguments correspond to the foil shapes, chord and span lengths, twist, dihedral and sweep angles.
# !!! info
#     A **wing section** consists of two foil profiles and their chord lengths and twist angles. Between them is their span length with associated _leading-edge_ dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for spans $b$, dihedrals $\delta$, sweeps $\Lambda$, and ``n+1`` entries for foils, chords $c$, and twists $\iota$, for some ``n \in \mathbb N``.
# 
# ![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)
airfoil = naca4((2,4,1,2))
wing = Wing(
    foils     = [ airfoil for i in 1:3 ],   # Foils
    chords    = [0.4, 0.2, 0.1],  # Chord lengths
    twists    = [0., 2., 5.],     # Twist angles (degrees)
    spans     = [1.0, 0.1],       # Section span lengths
    dihedrals = [0., 60.],        # Dihedral angles (degrees)
    sweeps    = [0., 30.],        # Sweep angles (degrees)
    chord_ratio   = 0.25,             # Sweep angle location w.r.t. 
                                  ## normalized chord lengths ∈ [0,1]
    symmetry  = true,             # Whether wing is symmetric
    ## flip      = false           # Whether wing is flipped in x-z plane
)

# The `symmetry` Boolean argument specifies whether the geometry should be reflected in the ``x``-``z`` plane. The `flip` Boolean argument specifies whether the coordinates should be flipped in the ``x``-``z`` plane. The `chord_ratio` argument specifies the chordwise-ratio of the sweep angles, e.g. 0. = leading edge sweep angle (default), 1. = trailing edge, 0.25 = quarter-chord.

# The following "getter" functions provide quantities of interest such as chord lengths, spans, twist, dihedral, and sweep angles.

# ```@repl howto
# foils(wing)
# chords(wing)
# spans(wing)
# twists(wing)
# dihedrals(wing)
# sweeps(wing) # Leading-edge sweep angles
# sweeps(wing, 0.25) # Quarter-chord sweep angles
# sweeps(wing, 1.0) # Trailing-edge sweep angles
# ```

# You can evaluate commonly defined properties such as the aspect ratio, projected area, etc., with the following functions.

# ```@repl howto
# aspect_ratio(wing)
# projected_area(wing)
# span(wing)
# taper_ratio(wing)
# mean_aerodynamic_chord(wing)
# ```

# The mean aerodynamic center calculation merits discussion. Its calculation procedure is depicted in the following diagram:
# 
# ![](https://raw.githubusercontent.com/HKUST-OCTAD-LAB/MECH3620Materials/main/pics/WingParams.svg)
mac25_w = mean_aerodynamic_center(wing) # 25% MAC by default
mac40_w = mean_aerodynamic_center(wing, 0.4) # at 40% MAC

# You can plot the wing with `Plots.jl` quite simply.
plot(wing, 
    0.25; # Aerodynamic center location in terms of chord ratio plot 
    mac = true, # Optional named argument to disable aerodynamic center plot
    zlim = (-0.5, 0.5) .* span(wing), 
    aspect_ratio = 1, 
    label = "Wing",
)

# At 40% MAC
plot!(wing, 
    0.40; # 40% MAC
    label = "40% MAC"
)

# There is also a convenient function for pretty-printing information (using [PrettyTables.jl](https://ronisbr.github.io/PrettyTables.jl/stable/)), in which the first argument takes the `Wing` type and the second takes a name (as a `String` or `Symbol`).
print_info(wing, "Wing")

#md # !!! tip
#md #     You can use [Accessors.jl](https://github.com/JuliaObjects/Accessors) to conveniently copy and modify properties of an existing object.

## Import Accessors
using Accessors

## Set only chords with other properties remaining identical.
new_wing = @set wing.chords = [0.4, 0.1, 0.05]

# ## Vortex Lattice Aerodynamic Analyses
# 
# How to run a generic aerodynamic analysis on a conventional aircraft configuration.
#
# ### Geometry 
# First we define the lifting surfaces. These can be a combination of `Wing` types constructed using the various methods available.
# !!! info
#     Support for fuselages and control surfaces is in progress.

# Horizontal tail
htail = WingSection(
    area       = 2.56 / 2,
    aspect     = 6.25,
    dihedral   = 0.0,
    sweep      = 15.0,
    taper      = 0.6,
    root_twist = 0.0,
    tip_twist  = 0.0,
    root_foil  = naca4(0,0,1,2),
    tip_foil   = naca4(0,0,0,9),
    position   = [5., 0., -0.1],
    angle      = 0.,
    axis       = [0., 1., 0.],
    symmetry   = true
);

# Vertical tail
vtail = WingSection(
    area       = 0.512 / 2,
    aspect     = 1.25,
    dihedral   = 0.0,
    sweep      = 8.0,
    taper      = 0.6,
    root_twist = 0.0,
    tip_twist  = 0.0,
    root_foil  = naca4(0,0,0,9),
    tip_foil   = naca4(0,0,0,9),
    position   = [5., 0., 0.],
    angle      = 90.,
    axis       = [1., 0., 0.]
)

# ### Meshing

# The `WingMesh` type takes a `Wing` type with a vector of integers consisting of the spanwise panel distribution corresponding to the number of sections, and an integer for the chordwise distribution.

## Wing meshes
wing_mesh = WingMesh(wing, [20,8], 10)

# !!! tip 
#     Optionally ,the type of spanwise spacing can be specified by the keyword `span_spacing` and providing types `Sine(), Cosine(), Uniform()` or a vector of the combination with length corresponding to the number of sections. For example, if you have two spanwise sections and the wing is symmetric, then the total number of spanwise sections is four. So the spanwise spacing argument would be, for example, `spanwise_spacing = [Uniform(), Cosine(), Cosine(), Uniform()]`

## Horizontal tail
htail_mesh = WingMesh(htail, [10], 5;
                      span_spacing = Cosine()) # Uniform(), Sine()

## Vertical tail
vtail_mesh = WingMesh(vtail, [10], 5); # Vertical tail

# ### Inviscid Analysis

# The inviscid 3D analysis uses a vortex lattice method. The `WingMesh` type allows you to generate horseshoe elements easily.
wing_horsies = make_horseshoes(wing_mesh)

# !!! note
#     You can also generate vortex ring elements by calling `make_vortex rings`, but the aerodynamic analysis is slightly slower.

# For multiple lifting surfaces, it is most convenient to define a single vector consisting of all the components' horseshoes using [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl).
aircraft = ComponentVector(
    wing = wing_horsies,
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

# To define reference values, use the following `References` type.

## Define reference values
refs = References(
    speed     = 150.0,
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing), 
    chord     = mean_aerodynamic_chord(wing), 
    location  = mean_aerodynamic_center(wing)
)

# The vortex lattice analysis can be executed with the vortex elements, freestream conditions, and reference values defined. An optional named argument is provided for enabling compressibility corrections at ``M > 0.3`` via the Prandtl-Glauert transformation.

## Solve system
system = solve_case(
    aircraft, fs, refs;
    compressible = true, # Compressibility corrections via Prandtl-Glauert transformation
    print = true, # Prints the results for only the aircraft
    print_components = true, # Prints the results for all components
)

# If needed, you can access the relevant component influence matrix values and boundary conditions with the following attributes, e.g.
system.influence_matrix[:wing] # Or :htail, :vtail
system.influence_matrix[:wing, :htail];

system.boundary_vector[:wing]

# ### Forces and Moments
# 
# The analysis computes the aerodynamic force coefficients $C_{\mathbf{F}}$ by surface pressure integration for nearfield forces. The aerodynamic moment coefficients $C_{\mathbf M}$ are calculated by computing the cross product of the moment arm of each control point from the reference point $\mathbf r_i - \mathbf r_{ref}$, i.e. $(C_{\mathbf{M}})_i = (\mathbf r_i - \mathbf r_{ref}) \times (C_{\mathbf{F}})_i$

# !!! note
#     You can specify the axis system for the nearfield forces, with choices of geometry, body, wind, or stability axes. The wind axes are used by default.

## Compute dynamics
ax = Wind() # Body(), Stability(), Geometry()

## Compute non-dimensional coefficients
CFs, CMs = surface_coefficients(system; axes = ax)

# You can access the corresponding values of the components' by the name provided in the `ComponentVector`.
CFs.wing

# Functions are available for directly retrieving the dimensionalized forces and moments.
Fs, Ms = surface_dynamics(system; axes = ax)

# A Trefftz plane integration is performed to compute farfield forces.
# 
# !!! note
#     The farfield forces are usually more accurate compared to nearfield forces.
# 
# To obtain the nearfield coefficients of the components (in wind axes by definition):
nfs = nearfield_coefficients(system)

# Similarly for the farfield coefficients of the components. 
ffs = farfield_coefficients(system)

# You can access the values corresponding to the components by the name used in the `ComponentArray` construction.
@show (nfs.wing, ffs.wing)

# To obtain the total nearfield and farfield force coefficients:
nf, ff = nearfield(system), farfield(system)

# You can also print the relevant information as a pretty table, if necessary.
print_coefficients(nfs.wing, ffs.wing, :wing)
print_coefficients(nf, ff, :aircraft)

# ## Aerodynamic Stability Analyses
# The derivatives of the aerodynamic coefficients with respect to the freestream values are obtained by automatic differentiation enabled by [ForwardDiff](https://github.com/JuliaDiff/ForwardDiff.jl). The following function evaluates the derivatives of the aerodynamic coefficients with respect to the freestream values. You can also optionally provide the axes for the reference frame of the coefficients.
dv_data = freestream_derivatives(
    system,
    print = true,  # Prints the results for only the aircraft
    print_components = true,  # Prints the results for all components
);

#md # !!! note
#md #     For efficiency, instead of calling `solve_case` to compute the forces and then computing the derivatives, you can directly call:
#md #     ```julia
#md #     dvs = freestream_derivatives(
#md #               aircraft, fs, refs,
#md #               compressible     = true, # Compressibility option
#md #               print            = true, # Prints the results for only the aircraft
#md #               print_components = true, # Prints the results for all components
#md #           )
#md #     ```