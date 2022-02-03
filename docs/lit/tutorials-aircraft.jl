# ## Objectives
#
# Here we will show you how to perform an aerodynamic analysis of a conventional aircraft. Specifically we will:
# 1. Define the geometries of a wing, horizontal tail and vertical tail.
# 2. Mesh and plot the geometries for numerical analyses.
# 3. Perform an aerodynamic analysis of this aircraft configuration at given freestream conditions and reference values.
# 4. Evaluate its drag polar for a given range of angles of attack.
# 5. Plot the spanwise loading distribution for the wing.

# For this, we will need to import some packages which will be convenient for plotting.
using AeroMDAO      # Main package
using Plots         # Plotting library
gr(dpi = 300)       # Plotting backend
using LaTeXStrings  # For LaTeX printing in plots

# ## Your First Wing
# 
# Here you will learn how to define a wing using an intuitive parametrization scheme. First, we define a `Vector` of `Foil`s.

## Define one airfoil
airfoil_1 = naca4((4,4,1,2))

## Define vector of airfoils
airfoils  = [ airfoil_1, naca4((0,0,1,2)) ]

#md # !!! note
#md #     Refer to the [Airfoil Aerodynamic Analysis](tutorials-airfoil.md) tutorial for an introduction to the `Foil` type.

# ### Parametrization
# The following function defines a symmetric wing and prints the relevant information. Named arguments corresponding to the foil shapes, chord and span lengths, twist, dihedral and leading-edge sweep angles. The following parametrization is used for the wing, presented for visual understanding.
# ![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)
wing = Wing(foils     = airfoils,    # Foil profiles
            chords    = [1.0, 0.6],  # Chord lengths
            twists    = [2.0, 0.0],  # Twist angles (degrees)
            spans     = [4.0],       # Section span lengths
            dihedrals = [5.],        # Dihedral angles (degrees)
            sweeps      = [5.])        # Leading-edge sweep angles (degrees)

#md # !!! info
#md #     See the [how-to guide](howto.md) on how to define an asymmetric wing.

# ### Visualization
# 
# The following function generates the coordinates of the wing's outline.
wing_outline = plot_wing(wing)

# Let's plot the geometry!
plt = plot(
           wing_outline[:,1], wing_outline[:,2], wing_outline[:,3], 
           label = "Wing",
           xaxis = "x", yaxis = "y", zaxis = "z",
           aspect_ratio = 1, 
           camera = (30, 45),
           zlim = (-0.1, span(wing) / 2),
          )

# ## Your First Vortex Lattice Analysis
# Now we would like to analyze the aerodynamics of this wing in conjunction with other lifting surfaces pertinent to aircraft.
# 
# ### Geometry
# We define the horizontal tail similarly to the wing. However, we also add additional position (by specifying a vector) and orientation attributes (by specifying an angle and axis of rotation) to place it at the desired location.

## Horizontal tail
htail = Wing(foils     = fill(naca4(0,0,1,2), 2),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweeps      = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# For the vertical tail, we simply replace `Wing` with `HalfWing` to define its shape.
## Vertical tail
vtail = HalfWing(foils     = fill(naca4(0,0,0,9), 2),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweeps      = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Let's visualize the geometry of the aircraft's configuration.
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

# ### Meshing and Assembly
# To perform the aerodynamic analysis, we will need to discretize our geometry into a _mesh_. The following `WingMesh` function constructs a mesh for you by providing a `HalfWing` or `Wing` type with specification of spanwise panels and chordwise panels. As the wing has only one spanwise section, we provide a vector with a single  integer entry for the spanwise panel distribution.
wing_mesh = WingMesh(wing, [12], 6) # (Wing, [Spanwise panels], Chordwise panels)

# Let's see what this discretization looks like on the camber distribution of the wing.

## Compute camber panel distribution
wing_cam_panels = camber_panels(wing_mesh)

## Generate plotting points
plt_wing_pans   = plot_panels(wing_cam_panels)

[ plot!(plt, panel, label = "", color = :lightblue) for panel in plt_wing_pans ]
plot!(plt)

# Similarly we define the meshes for the other surfaces and plot them.
htail_mesh = WingMesh(htail, [6], 4)
vtail_mesh = WingMesh(vtail, [4], 3)

[ plot!(plt, panel, label = "", color = :orange) 
    for panel in plot_panels(camber_panels(htail_mesh)) ]
[ plot!(plt, panel, label = "", color = :lightgreen) 
    for panel in plot_panels(camber_panels(vtail_mesh)) ]
plot!(plt)

# For the analysis, you have to assemble the meshes into a `ComponentArray/Vector`.
aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         )

# You can define the freestream condition as follows, by providing the angles of attack $\alpha$ and sideslip $\beta$ in degrees with a rotation vector $\Omega$.
fs  = Freestream(
                 alpha = 3.0, # degrees
                 beta  = 0.0, # degrees
                 omega = [0., 0., 0.]
                );

# You can define the reference values for the speed, area, span, chord, density, and location  as follows.
refs = References(speed    = 1.0, 
                  area     = projected_area(wing),
                  span     = span(wing),
                  chord    = mean_aerodynamic_chord(wing),
                  density  = 1.225,
                  location = mean_aerodynamic_center(wing));

# You can run the aerodynamic analysis by providing the aircraft configuration, freestream, and reference values. Optionally you can also print the results.
system = solve_case(
        aircraft, fs, refs;
        print            = true, # Prints the results for only the aircraft
        print_components = true, # Prints the results for all components
    )

# You can obtain the aerodynamic coefficients from this system. The nearfield aerodynamic force and moment coefficients are ordered as $(C_{D_i}, C_Y, C_L, C_\ell, C_m, C_n)$.
nf = nearfield(system)

# !!! tip 
#     Refer to the [how-to guide](howto.md) to see how to compute the aerodynamic coefficients of each component and perform stability analyses.

# ### Drag Polar

# Now let's analyze the drag polar of this aircraft configuration by varying the angle of attack and collecting the induced drag coefficient $C_{D_i}$.

## Define function to compute system varying with angle of attack.
vary_alpha(aircraft, α, refs) = solve_case(aircraft, Freestream(alpha = α), refs)

## Run loop
αs      = -5:0.5:5
systems = [ vary_alpha(aircraft, α, refs) for α in αs ]
## Cleaner: map(α -> vary_alpha(...), αs)

## Get coefficients
coeffs = nearfield.(systems)
CDis   = [ c[1] for c in coeffs ]
CLs    = [ c[3] for c in coeffs ];

# Let's plot the drag polar!
plot(CDis, CLs, 
     label  = "",
     xlabel = L"C_{D_i}",
     ylabel = L"C_L",
     title  = "Drag Polar",
     ls     = :solid)

# Let's also take a look at the variations of all the coefficients.

## Concatenate results into one array
data = permutedims(reduce(hcat, [α; c] for (α, c) in zip(αs, coeffs)))

## Plot
plot(data[:,1], round.(data[:,2:end], digits = 4), 
     layout = (3,2),
     xlabel = L"\alpha",
     ylabel = [L"C_{D_i}" L"C_Y" L"C_L" L"C_\ell" L"C_m" L"C_n" ],
     labels = "",
     size   = (800, 600)
    )

# > **Tip:** You can convert this into a DataFrame for convenient reference.
# > ```julia
# > using DataFrames, StatsPlots
# > data = DataFrame([ xs for xs in zip(αs, coeffs) ])
# > rename!(data, [:α, :CD, :CY, :CL, :Cl, :Cm, :Cn])
# > @df data plot()
# > ```

# ### Spanwise Loading

# You can compute the aerodynamic coefficients on the panels from the system.
CFs, CMs = surface_coefficients(system)

## Get panels along the chord-lines of the wing.
wing_panels = chord_panels(wing_mesh)

## Compute spanwise loads
span_loads  = spanwise_loading(wing_panels, CFs.wing, projected_area(wing))
CL_loads    = vec(sum(system.circulations.wing, dims = 1)) / (0.5 * refs.speed * refs.chord);

## Plot spanwise loadings
plot_CD = plot(span_loads[:,1], span_loads[:,2], label = :none, ylabel = "CDi")
plot_CY = plot(span_loads[:,1], span_loads[:,3], label = :none, ylabel = "CY")
plot_CL = begin
            plot(span_loads[:,1], span_loads[:,4], label = :none, xlabel = "y", ylabel = "CL")
            plot!(span_loads[:,1], CL_loads, label = "Normalized", xlabel = "y")
          end
plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))
