# Aerodynamics

```@setup aeromdao
using AeroMDAO
using Plots
using LaTeXStrings
gr()
```

The following sections describe different use-cases of AeroMDAO.

## Doublet-Source Panel Method

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
airfoil = naca4((digits; sharp_trailing_edge = true))  #hide
V, α = 1.0, 5.0             #hide
uniform = Uniform2D(V, α)   #hide
cl, cls, cms, cps, panels = solve_case(airfoil, 
                                       uniform;
                                       viscous = false,
                                       sources = false, 
                                       wake_length = 1e3,
                                       wake_panels = 100,
                                       num_panels = 80)

upper, lower = split_surface(airfoil); #hide
x_upper, y_upper = upper[:,1], upper[:,2] #hide
x_lower, y_lower = lower[:,1], lower[:,2] #hide

println("Cl: $cl")
println("Σᵢ Clᵢ: $(sum(cls))")
println("Σᵢ Cmᵢ: $(sum(cms))")
```

```@example aeromdao
digits  = (0,0,1,2); #hide
airfoil = naca4((digits; sharp_trailing_edge = true))  #hide
V, α = 1.0, 5.0             #hide
uniform = Uniform2D(V, α)   #hide
cl, cls, cms, cps, panels = solve_case(airfoil, #hide
                                       uniform; #hide
                                       viscous = false, #hide
                                       sources = false,  #hide
                                       wake_length = 1e3, #hide
                                       wake_panels = 100, #hide
                                       num_panels = 80) #hide

upper, lower = split_surface(airfoil); #hide
x_upper, y_upper = upper[:,1], upper[:,2] #hide
x_lower, y_lower = lower[:,1], lower[:,2] #hide

get_surface_values(panels, vals, surf = "lower") = [ (collocation_point(panel)[1], val) for (val, panel) in zip(vals, panels) if panel_location(panel) == surf ] #hide
cp_lower = get_surface_values(panels, cps, "lower") #hide
cp_upper = get_surface_values(panels, cps, "upper") #hide
plot(yflip = true, xlabel=L"(x/c)", ylabel = L"$C_p$") #hide
plot!(cp_upper, ls = :dash, lw = 2, c = :cornflowerblue, label = "Upper") #hide
plot!(cp_lower, ls = :dash, lw = 2, c = :orange, label = "Lower") #hide
plot!(x_upper, -y_upper,  #hide
        ls = :solid, lw = 2, c = :cornflowerblue, label = "NACA $(digits...) Upper") #hide
plot!(x_lower, -y_lower,  #hide
        ls = :solid, lw = 2, c = :orange, label = "NACA $(digits...) Lower") #hide
```

## Vortex Lattice Method

The vortex lattice method used in AeroMDAO follows Mark Drela's *Flight Vehicle Aerodynamics*.

### Meshing

The geometry "engine" generates panels for horseshoes and the camber distribution using the airfoil data in the definition of the wing.

```julia
WingMesh(wing                    :: AbstractWing,     # Wing type
         span_num                :: Vector{Integer},  # Number of spanwise panels
         chord_num               :: Integer,          # Number of chordwise panels
         span_spacing = Cosine() :: AbstractSpacing   # Spanwise distribution
        )
```

```julia
wing_mesh = WingMesh(wing, [12, 3], 6);
wing_mesh.camber_mesh
```

### Aircraft Definition

Aircraft analysis by definition of multiple lifting surfaces using the `HalfWing` or `Wing` types are also supported, but is slightly more complex due to the increases in user specifications. Particularly, you will have to mesh the different components by yourself by specifying the number of chordwise panels, spanwise panels and their associated spacings.

This method returns the horseshoe panels for the analysis, and the associated normal vectors based on the camber distribution of the wing. Consider a "vanilla" aircraft, in which you have a wing, horizontal tail, and vertical tail.

```julia
# Wing
wing = Wing(foils     = Foil.(fill(naca4(2,4,1,2), 2)),
            chords    = [1.0, 0.6],
            twists    = [2.0, 0.0],
            spans     = [4.0],
            dihedrals = [5.],
            sweeps      = [5.]);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(foils     = Foil.(fill(naca4(0,0,1,2), 2)),
             chords    = [0.7, 0.42],
             twists    = [0.0, 0.0],
             spans     = [1.25],
             dihedrals = [0.],
             sweeps      = [6.39],
             position  = [4., 0, 0],
             angle     = -2.,
             axis      = [0., 1., 0.])

# Vertical tail
vtail = HalfWing(foils     = Foil.(fill(naca4(0,0,0,9), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweeps      = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Print info
print_info(wing,  "Wing")
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

The horseshoes for the analysis are generated using the `make_horseshoes` function.

```julia
wing_horsies = make_horseshoes(wing_mesh),
```

To run analyses with multiple components, you will have to assemble this information using the exported `ComponentVector` (built by [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)). This names the components using a symbolic tag, and all the components are treated as a single vector.

```julia
aircraft = ComponentVector(
                           wing  = wing_horsies,
                           htail = make_horseshoes(htail_mesh),
                           vtail = make_horseshoes(vtail_mesh)
                          );
```

### Freestream Condition

This geometry is analysed at given freestream angles of attack and sideslip (in degrees), which define the boundary conditions of the problem.

```julia
fs  = Freestream(alpha = 0.0, 
                 beta  = 0.0, 
                 omega = [0., 0., 0.]);
```

### References

For an analysis, you require ambient reference conditions. Specifically, you need the density $\rho$ and a reference location $r_\text{ref}$ for calculating forces and moments. These dynamics are non-dimensionalized using the wing's characteristics as reference.

```julia
refs = References(speed    = 1.0, 
                  density  = 1.225,
                  area     = projected_area(wing),
                  span     = span(wing),
                  chord    = mean_aerodynamic_chord(wing),
                  location = mean_aerodynamic_center(wing))
```

### Analysis

```julia
system = solve_case(
                    aircraft   :: Vector{Horseshoe}, # Horseshoes on the aircraft surface
                    fs         :: Freestream,        # Freestream values
                    refs       :: References;        # Reference values
                    print            = true,         # Prints the results for only the aircraft
                    print_components = true,         # Prints the results for all components
                   );
```

It returns a system which can be used to determine quantities of interest such as the dynamics.

### Dynamics

The analysis computes the aerodynamic forces by surface pressure integration for nearfield forces. You can specify the axis system for the nearfield forces, with choices of `Geometry(), Body(), Wind(), Stability()`. The wind axes are used by default.

```julia
# Compute dynamics
ax       = Geometry() # Body(), Stability(), Wind()
CFs, CMs = surface_coefficients(system; axes = ax)
Fs       = surface_forces(system)
Fs, Ms   = surface_dynamics(system; axes = ax) 
```

A Trefftz plane integration is employed to obtain farfield forces.

!!! note
    The farfield forces are usually more accurate compared to nearfield forces, as the components do not interact as in the evaluation of the Biot-Savart integrals for the latter.

To obtain the nearfield and farfield coefficients of the components (in wind axes by definition):

```julia
nfs = nearfield_coefficients(system)
ffs = farfield_coefficients(system)
```

To obtain the total nearfield and farfield force coefficients:

```julia
nf  = nearfield(system) 
ff  = farfield(system)
```

### Derivatives

AeroMDAO uses the [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) package, which leverages forward-mode automatic differentiation to obtain derivatives of certain quantities of interest. Particularly, with respect to the angles of attack and sideslip, and the non-dimensionalized roll rates at low computational cost for stability analyses. A unique feature which is not easily obtained from other implementations is the stability derivatives of each component. To obtain these quantities, simply replace `solve_case()` with `solve_case_derivatives()`, which will return data consisting of the nearfield, farfield and stability derivative coefficients.

!!! note
    This is due to limitations of using closures in `ForwardDiff`, and hence only these quantities are provided using this function.

```julia
solve_case_derivatives(
                       aircraft, fs, refs;
                       print            = true,    # Prints the results for only the aircraft
                       print_components = true,    # Prints the results for all components
                      );
```

You can pretty-print the stability derivatives with the following function, whose first argument again provides the name of the wing:

```julia
print_derivatives(dv_data.wing, "My Wing")
```