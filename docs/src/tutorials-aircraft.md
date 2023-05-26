```@meta
EditURL = "<unknown>/docs/lit/tutorials-aircraft.jl"
```

## Objectives

Here we will show you how to perform an aerodynamic analysis of a conventional aircraft.
> **Recipe**
> 1. Define the geometries of a wing, horizontal tail and vertical tail.
> 2. Mesh and plot the geometries for numerical analyses.
> 3. Perform an aerodynamic analysis of this aircraft configuration at given freestream conditions and reference values.
> 4. Evaluate its drag polar for a given range of angles of attack.
> 5. Plot the spanwise loading distribution for the wing.

For this, we will need to import some packages which will be convenient for plotting.

````@example tutorials-aircraft
using AeroFuse          # Main package
using Plots             # Plotting library
gr(                     # Plotting backend
    size = (800,600),   # Size
    dpi = 300,          # Resolution
    palette = :Dark2_8  # Color scheme
)
using LaTeXStrings      # For LaTeX printing in plots
````

## Your First Wing

Here you will learn how to define a wing using an intuitive parametrization scheme. First, we define a `Vector` of `Foil`s.

````@example tutorials-aircraft
# Define one airfoil
airfoil_1 = naca4(4,4,1,2)

# Define vector of airfoils
airfoils  = [ airfoil_1, naca4((0,0,1,2)) ]
````

!!! note
    Refer to the [Airfoil Aerodynamic Analysis](tutorials-airfoil.md) tutorial for an introduction to the `Foil` type.

### Parametrization
The following function defines a symmetric, single-section wing and prints the relevant information.
!!! note
    All relevant units are presented in metric.

````@example tutorials-aircraft
wing = Wing(
    foils     = airfoils,    # Foil profiles
    chords    = [1.0, 0.6],  # Chord lengths
    twists    = [2.0, 0.0],  # Twist angles (degrees)
    spans     = [4.0],       # Section span lengths
    dihedrals = [5.],        # Dihedral angles (degrees)
    sweeps    = [5.],        # Sweep angles (degrees)
    chord_ratio   = 0.,          # Sweep angle location w.r.t.
                             # normalized chord lengths ∈ [0,1]
    symmetry  = true,        # Whether wing is symmetric
    # flip      = false      # Whether wing is reflected
)
````

In this case, the root chord is $1.0~m$ and the tip chord is $0.6~m$, and so on.

!!! info
    See the [How-to Guide](howto.md) for the various ways of constructing wings.

### Visualization

You can plot the lifting surface by calling `plot` with `Plots.jl`.

````@example tutorials-aircraft
plt = plot(
    xlabel = L"x", ylabel = L"y", zlabel = L"z",
    aspect_ratio = 1,
    camera = (30, 30),
    zlim = (-0.5, 0.5) .* span(wing),
)

plot!(plt, wing, label = "Wing")
````

## Your First Vortex Lattice Analysis
Now we would like to analyze the aerodynamics of this wing in conjunction with other lifting surfaces pertinent to aircraft.

### Geometry
We define the horizontal tail similarly to the wing. However, we also add additional position (by specifying a vector) and orientation attributes (by specifying an angle and axis of rotation) to place it at the desired location.

````@example tutorials-aircraft
# Horizontal tail
htail = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25],
    dihedrals = [0.],
    sweeps    = [6.39],
    chord_ratio   = 0.,
    position  = [4., 0, 0],
    angle     = -2.,
    axis      = [0., 1., 0.],
    symmetry  = true
)
````

For the vertical tail, we simply set `symmetry = false` to define its shape.

````@example tutorials-aircraft
# Vertical tail
vtail = Wing(
    foils     = fill(naca4(0,0,0,9), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.0],
    dihedrals = [0.],
    sweeps    = [7.97],
    chord_ratio   = 0.,
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.]
)
````

Let's visualize the geometry of the aircraft's configuration.

````@example tutorials-aircraft
plot!(plt, htail, label = "Horizontal Tail")
plot!(plt, vtail, label = "Vertical Tail")
````

### Meshing
To perform the aerodynamic analysis, we will need to discretize our geometry into a _mesh_. The following `WingMesh` function constructs a mesh for you by providing a `HalfWing` or `Wing` type with specification of spanwise panels and chordwise panels. As the wing has only one spanwise section, we provide a vector with a single  integer entry for the spanwise panel distribution.

````@example tutorials-aircraft
wing_mesh = WingMesh(wing, [20], 12) # (Wing, [Spanwise panels], Chordwise panels)
````

Let's see what this discretization looks like on the camber distribution of the wing.

````@example tutorials-aircraft
# Compute camber panel distribution
plot!(plt, wing_mesh, label = "")
````

Similarly we define the meshes for the other surfaces and plot them.

````@example tutorials-aircraft
htail_mesh = WingMesh(htail, [8], 6)
vtail_mesh = WingMesh(vtail, [6], 4)

plot!(plt, htail_mesh, label = "")
plot!(plt, vtail_mesh, label = "")
````

### Aerodynamic Analysis
For the analysis, we need to generate a `Horseshoe` type, corresponding to horseshoe singularity elements used in the vortex lattice method. This is generated as follows:

````@example tutorials-aircraft
make_horseshoes(wing_mesh)
````

To perform the aerodynamic analysis, you have to assemble these horseshoes for each surface into a `ComponentVector`.

````@example tutorials-aircraft
aircraft = ComponentVector(
    wing  = make_horseshoes(wing_mesh),
    htail = make_horseshoes(htail_mesh),
    vtail = make_horseshoes(vtail_mesh)
)
````

You can define the freestream condition as follows, by providing the angles of attack $\alpha$ and sideslip $\beta$ in degrees with a rotation vector $\Omega$.

````@example tutorials-aircraft
fs  = Freestream(
    alpha = 3.0, # degrees
    beta  = 0.0, # degrees
    omega = [0., 0., 0.]
);
nothing #hide
````

You can define the reference values for the speed, area, span, chord, density, and location  as follows.

````@example tutorials-aircraft
refs = References(
    speed    = 1.0,
    area     = projected_area(wing),
    span     = span(wing),
    chord    = mean_aerodynamic_chord(wing),
    density  = 1.225,
    location = mean_aerodynamic_center(wing)
);
nothing #hide
````

You can run the aerodynamic analysis by providing the aircraft configuration, freestream, and reference values. Optionally you can also print the results.

````@example tutorials-aircraft
system = solve_case(
    aircraft, fs, refs;
    compressible     = false, # Compressibility option
    # print            = true, # Prints the results for only the aircraft
    # print_components = true, # Prints the results for all components
)
````

You can obtain the aerodynamic coefficients from this system. The nearfield aerodynamic force and moment coefficients are ordered as $(C_{D_i}, C_Y, C_L, C_\ell, C_m, C_n)$.

````@example tutorials-aircraft
nf = nearfield(system)
````

The force coefficients are printed as $(C_X, C_Y, C_Z)$ for general axis systems; wind axes are used in the `nearfield` function.

!!! tip
    Refer to the [how-to guide](howto.md) to see how to compute the aerodynamic coefficients of each component and perform stability analyses.

A convenience method is also provided for plotting streamlines from the leading edge of each surface.

````@example tutorials-aircraft
plot!(plt,
    system,     # VortexLattice System
    wing_mesh,  # Lifting surface (or mesh)
    span = 4,   # Number of streamlines per spanwise section
    dist = 10,  # Distance of streamlines (m)
    # lc = :green, # Color of streamlines
)

plot!(plt, system, htail_mesh, span = 3, lc = :cyan) # For horizontal tail
plot!(plt, system, vtail, span = 2, lc = :cyan) # For vertical tail
````

### Drag Polar

Now let's analyze the drag polar of this aircraft configuration by varying the angle of attack and collecting the induced drag coefficient $C_{D_i}$.

````@example tutorials-aircraft
# Define function to compute system varying with angle of attack.
vary_alpha(aircraft, α, refs) = VortexLatticeSystem(aircraft, Freestream(alpha = α), refs)

# Run loop
αs      = -5:0.5:5
systems = [ vary_alpha(aircraft, α, refs) for α in αs ]
# Cleaner: map(α -> vary_alpha(...), αs)

# Get coefficients
coeffs = nearfield.(systems)
CDis   = [ c[1] for c in coeffs ]
CLs    = [ c[3] for c in coeffs ];
nothing #hide
````

Let's plot the drag polar!

````@example tutorials-aircraft
plot(CDis, CLs,
     label  = "",
     xlabel = L"C_{D_i}",
     ylabel = L"C_L",
     title  = "Drag Polar",
     ls     = :solid)
````

Let's also take a look at the variations of all the coefficients.

````@example tutorials-aircraft
# Concatenate results into one array
data = permutedims(
    mapreduce(hcat, systems) do sys
        [sys.freestream.alpha; nearfield(sys) ]
    end
)

# Plot
plot(
    rad2deg.(data[:,1]),  # Angle of attack
    round.(data[:,2:end], digits = 4), # Aerodynamic coefficients
    layout = (3,2),
    xlabel = L"\alpha",
    ylabel = [L"C_{D_i}" L"C_Y" L"C_L" L"C_\ell" L"C_m" L"C_n"],
    labels = "",
)
````

> **Tip:** You can convert this into a DataFrame for convenient reference.

````@example tutorials-aircraft
using DataFrames
df = DataFrame(round.(data, digits = 6), [:α, :CX, :CY, :CZ, :Cl, :Cm, :Cn])
````

### Spanwise Loading

You can compute the aerodynamic coefficients on the panels from the system.

````@example tutorials-aircraft
CFs, CMs = surface_coefficients(system)

# Compute spanwise loads
span_loads = spanwise_loading(wing_mesh, system.reference, CFs.wing, system.circulations.wing)

# Plot spanwise loadings
plot_CD = plot(span_loads[:,1], span_loads[:,2], label = :none, ylabel = L"C_{D_i}")
plot_CY = plot(span_loads[:,1], span_loads[:,3], label = :none, ylabel = L"C_Y")
plot_CL = begin
            plot(span_loads[:,1], span_loads[:,4], label = :none, xlabel = L"y", ylabel = L"C_L")
            plot!(span_loads[:,1], span_loads[:,5], label = L"Γ/ρVc", xlabel = L"y")
        end
plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

