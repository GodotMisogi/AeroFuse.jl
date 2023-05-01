# ## Objectives
#
# Here we will show you how to perform an aerodynamic stability analysis of a conventional aircraft. Here, we'll attempt to replicate the design of a Boeing 777. It won't be a realistic replication, but the overall geometry will match to a certain extent.
# 
# ![](https://images-wixmp-ed30a86b8c4ca887773594c2.wixmp.com/f/dc763bf2-302c-46be-8a52-4cb7c11598e5/d74vi3c-372cf93b-f4ad-4046-85e3-49f667d3c55a.png?token=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJ1cm46YXBwOjdlMGQxODg5ODIyNjQzNzNhNWYwZDQxNWVhMGQyNmUwIiwiaXNzIjoidXJuOmFwcDo3ZTBkMTg4OTgyMjY0MzczYTVmMGQ0MTVlYTBkMjZlMCIsIm9iaiI6W1t7InBhdGgiOiJcL2ZcL2RjNzYzYmYyLTMwMmMtNDZiZS04YTUyLTRjYjdjMTE1OThlNVwvZDc0dmkzYy0zNzJjZjkzYi1mNGFkLTQwNDYtODVlMy00OWY2NjdkM2M1NWEucG5nIn1dXSwiYXVkIjpbInVybjpzZXJ2aWNlOmZpbGUuZG93bmxvYWQiXX0.bS5c5rkhqB2yoaOmIeRut7TgVsqgnPIfMOBSgYOO-TI)
# **Source**: [boeingboeing2, DeviantArt](https://www.deviantart.com/boeingboeing2/art/Boeing-777-vector-431451480)
#md # !!! note
#md #     Refer to the [Aircraft Aerodynamic Analysis](tutorials-aircraft.md) tutorial before studying this tutorial.
# > **Recipe**
# > 1. Define the geometries of a fuselage, wing, horizontal tail and vertical tail.
# > 2. Perform an aerodynamic analysis of this aircraft configuration at given freestream conditions and reference values.
# > 4. Evaluate the derivatives of the aerodynamic coefficients with respect to the freestream conditions.
# > 5. Evaluate the quantities of interest for aerodynamic stability.

# Let's import the relevant packages.
using AeroFuse          # Main package
using Plots             # Plotting library
gr(                     # Plotting backend
    size = (800,600),   # Size
    dpi = 300,          # Resolution
    palette = :Dark2_8  # Color scheme
) 
using LaTeXStrings      # For LaTeX printing in plots

# ## Aircraft Geometry

# First, let's define the fuselage. Here we'll define it by combining a cylindrical definition for the cabin with [hyperelliptical](https://mathworld.wolfram.com/Hyperellipse.html) curves for the nose and rear.

## Fuselage definition
fuse = HyperEllipseFuselage(
    radius = 3.04,          # Radius, m
    length = 63.7,          # Length, m
    x_a    = 0.15,          # Start of cabin, ratio of length
    x_b    = 0.7,           # End of cabin, ratio of length
    c_nose = 2.0,           # Curvature of nose
    c_rear = 1.2,           # Curvature of rear
    d_nose = -0.5,          # "Droop" or "rise" of nose, m
    d_rear = 1.0,           # "Droop" or "rise" of rear, m
    position = [0.,0.,0.]   # Set nose at origin, m
)

## Compute geometric properties
ts = 0:0.1:1                # Distribution of sections for nose, cabin and rear
S_f = wetted_area(fuse, ts) # Surface area, m²
V_f = volume(fuse, ts)      # Volume, m³

## Get coordinates of rear end
fuse_end = fuse.affine.translation + [ fuse.length, 0., 0. ]

# Now, let's define the lifting surfaces. We'll download a supercritical airfoil for the wing section; note that this is not the same one as used in the Boeing 777-200LR. We'll also define a two-section wing.

## Define one airfoil
foil_w = read_foil(download("http://airfoiltools.com/airfoil/seligdatfile?airfoil=rae2822-il"))

## Define vector of airfoils
foils = [ foil_w, foil_w, naca4((0,0,1,2)) ]

## Wing
wing = Wing(
    foils       = foils,                        # Airfoils
    chords      = [14.0, 9.73, 1.43561],        # Chord lengths
    spans       = [14.0, 46.9] / 2,             # Span lengths
    dihedrals   = fill(6, 2),                   # Dihedral angles (deg)
    sweeps      = fill(35.6, 2),                # Sweep angles (deg )
    w_sweep     = 0.,                           # Leading-edge sweep
    position    = [19.51, 0., -2.5],            # Position
    symmetry    = true                          # Symmetry
)

b_w = span(wing)
S_w = projected_area(wing)
c_w = mean_aerodynamic_chord(wing)
x_w, y_w, z_w = mac_w = mean_aerodynamic_center(wing)

# For reference, let's plot what we have so far.
p1 = plot(
    xaxis = L"x", yaxis = L"y", zaxis = L"z",
    aspect_ratio = 1, 
    zlim = (-0.5, 0.5) .* span(wing),
    camera = (30,30)
)

plot!(fuse, label = "Fuselage", alpha = 0.6)
plot!(wing, label = "Wing")

# ### Stabilizers
# Now, let's add the stabilizers. First, the horizontal tail.
htail = WingSection(
    area        = 101,  # Area (m²)
    aspect      = 4.2,  # Aspect ratio
    taper       = 0.4,  # Taper ratio
    dihedral    = 7.,   # Dihedral angle (deg)
    sweep       = 35.,  # Sweep angle (deg)
    w_sweep     = 0.,   # Leading-edge sweep
    root_foil   = naca4(0,0,1,2),
    symmetry    = true,
    
    ## Orientation
    angle       = -2,           # Incidence angle (deg)
    axis        = [0., 1., 0.], # Axis of rotation, y-axis
    position    = [ fuse_end.x - 8., 0., 0.],
)

b_h = span(htail)
S_h = projected_area(htail)
c_h = mean_aerodynamic_chord(htail)
x_h, y_h, z_h = mac_h = mean_aerodynamic_center(htail)

# Now the vertical tail.
vtail = WingSection(
    area        = 56.1, # Area (m²)
    aspect      = 1.5,  # Aspect ratio
    taper       = 0.4,  # Taper ratio
    sweep       = 44.4, # Sweep angle (deg)
    w_sweep     = 0.,   # Leading-edge sweep
    root_foil   = naca4(0,0,0,9),
    
    ## Orientation
    angle       = 90.,       # To make it vertical
    axis        = [1, 0, 0], # Axis of rotation, x-axis
    position    = htail.affine.translation - [2.,0.,0.]
) # Not a symmetric surface

b_v = span(vtail)
S_v = projected_area(vtail)
c_v = mean_aerodynamic_chord(vtail)
x_v, y_v, z_v = mac_v = mean_aerodynamic_center(vtail)

# Let's mesh and plot the lifting surfaces.
wing_mesh = WingMesh(wing, [8,16], 10, 
    span_spacing = fill(Uniform(), 4) # Number of spacings = number of spanwise stations
    ## (including symmetry)
)

htail_mesh = WingMesh(htail, [10], 8)
vtail_mesh = WingMesh(vtail, [8], 6)

## Plot meshes
plt = plot(
    xaxis = L"x", yaxis = L"y", zaxis = L"z",
    aspect_ratio = 1, 
    zlim = (-0.5, 0.5) .* span(wing_mesh),
    camera = (30, 45),
)
plot!(fuse, label = "Fuselage", alpha = 0.6)
plot!(plt, wing_mesh, label = "Wing")
plot!(plt, htail_mesh, label = "Horizontal Tail")
plot!(plt, vtail_mesh, label = "Vertical Tail")


# ## Aerodynamic Analysis

# Now, let's generate the horseshoe system for the aircraft.
aircraft = ComponentVector(
    wing  = make_horseshoes(wing_mesh),
    htail = make_horseshoes(htail_mesh),
    vtail = make_horseshoes(vtail_mesh)
)

#md # !!! warning "Alert!"
#md #      Note that the fuselage is not included in the analysis for now, and is only included for plotting. Support for its aerodynamic analysis will be added soon.

# Now, let's define the freestream conditions.
fs  = Freestream(
    alpha = 3.0, # degrees
    beta  = 0.0, # degrees
    omega = [0., 0., 0.]
);

# Similarly, define the reference values. Here, the reference flight condition will be set to Mach number $M = 0.84$.
M = 0.84 # Mach number
refs = References(
    sound_speed = 330.,
    speed    = M * 330., 
    density  = 1.225,
    span     = b_w,
    area     = S_w,
    chord    = c_w,
    location = mac_w
);

# Let's run the aerodynamic analysis first.
system = solve_case(
    aircraft, fs, refs;
    compressible     = true, # Compressibility option
    ## print            = true, # Prints the results for only the aircraft
    ## print_components = true, # Prints the results for all components
)

# !!! info 
#     You may receive a warning that the results are incorrect; this is due to the limitation of the vortex lattice method being able to primarily analyze only subsonic flows under the physical assumptions.

# The freestream derivatives can be obtained by passing the resultant system as follows:
dvs = freestream_derivatives(
    system,                     # VortexLatticeSystem
    axes             = Wind(),  # Specify axis system for nearfield forces (wind by default)
    ## print            = true,    # Prints the results for only the aircraft
    ## print_components = true,    # Prints the results for all components
    ## farfield         = true,    # Print farfield derivatives
);

# You can access the derivatives of each lifting surface based on the keys defined in the `ComponentVector`.
ac_dvs = dvs.aircraft

# These quantities are the force and moment coefficients $(C_X, C_Y, C_Z, C_l, C_m, C_n, C_{D_{i,ff}}, C_{Y_{ff}} C_{L_{ff}})$ generated from the nearfield and farfield analyses, and their derivatives respect to the Mach number $M$, freestream angles of attack and sideslip $(\alpha, \beta)$, and the non-dimensional angular velocity rates in stability axes $(\bar{p}, \bar{q}, \bar{r})$. The keys corresponding to the freestream derivatives should be evident:
keys(dvs.aircraft)

# These can be accessed either like a dictionary, or by 'dot' syntax.
ac_dvs[:CZ_al], ac_dvs.CZ_al, ac_dvs.CL_al # Lift coefficient derivative wrt. alpha 

# Note that the nearfield forces and moments $(C_X, C_Y, C_Z, C_l, C_m, C_n)$ depend on the axis system used ($C_Z$ is not lift if body axes are used!). You can also pretty-print the derivatives for each surface.
print_derivatives(dvs.aircraft, "Aircraft", farfield = true)
print_derivatives(dvs.wing, "Wing", farfield = true)
print_derivatives(dvs.htail, "Horizontal Tail", farfield = true)
print_derivatives(dvs.vtail, "Vertical Tail", farfield = true)

# ## Static Stability Analysis

# You can evaluate the static stability of the aircraft using the various quantities computed in this process, following standard computations and definitions in aircraft design and stability analysis.
l_h = x_h - x_w                 # Horizontal tail moment arm
V_h = S_h / S_w * l_h / c_w     # Horizontal tail volume coefficient

l_v = x_v - x_w                 # Vertical tail moment arm
V_v = S_v / S_w * l_v / b_w     # Vertical tail volume coefficient

x_cp = -refs.chord * ac_dvs.Cm / ac_dvs.CZ # Center of pressure

x_np_lon = -refs.chord * ac_dvs.Cm_al / ac_dvs.CZ_al # Neutral point

## Position vectors
x_np, y_np, z_np = r_np = refs.location + [ x_np_lon; zeros(2) ]  # Neutral point
x_cp, y_cp, z_cp = r_cp = refs.location + [ x_cp; zeros(2) ]  # Center of pressure

@info "Horizontal TVC                V_h:" V_h
@info "Vertical TVC                  V_v:" V_v
@info "Wing Aerodynamic Center  x_ac (m):" x_w
@info "Neutral Point            x_np (m):" x_np
@info "Center of Pressure       x_cp (m):" x_cp

# Let's plot the relevant quantities.
stab_plt = plot(
    xaxis = L"x", yaxis = L"y", zaxis = L"z",
    aspect_ratio = 1, 
    zlim = (-0.5, 0.5) .* span(wing_mesh),
    camera = (30,60),
)
plot!(fuse, label = "Fuselage", alpha = 0.6)
plot!(stab_plt, wing_mesh, label = "Wing", mac = false)
plot!(stab_plt, htail_mesh, label = "Horizontal Tail", mac = false)
plot!(stab_plt, vtail_mesh, label = "Vertical Tail", mac = false)

scatter!(Tuple(r_np), color = :orange, label = "Neutral Point")
scatter!(Tuple(r_cp), color = :brown, label = "Center of Pressure")