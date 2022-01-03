#  Guide

```@setup aeromdao
using AeroMDAO
using StaticArrays
using Plots
gr()
```

AeroMDAO provides convenient methods for developing studies in aerodynamics.

---

## Airfoil Parametrization

AeroMDAO provides some basic parametrizations commonly used for airfoils, with conversion to the camber-thickness representation and vice versa.

*NACA 4-digit Airfoil Parametrization*

```julia
naca4(digits :: NTuple{4, Real})
```

```julia
airfoil = naca4((2,4,1,2))
```

*Kulfan Class Shape Transformation (CST) Method*
```julia
kulfan_CST(alpha_u      :: Vector{Real},    # Upper surface parameters
           alpha_l      :: Vector{Real},    # Lower surface parameters
           dzs          :: NTuple{2, Real}, # Upper and lower trailing edge points
           coeff_LE = 0 :: Real,            # Leading-edge modification coefficient
           n = 40       :: Integer)         # Number of points on each surface
```

```julia
alpha_u = [0.1, 0.3, 0.2, 0.15, 0.2]
alpha_l = [-0.1, -0.1, -0.1, -0.001, -0.02]
dzs = (1e-4, 1e-4)
foil = kulfan_CST(alpha_u, alpha_l, dzs, 0.2)
```

*Camber-Thickness Representation*
```julia
foil_camthick(coords    :: Array{2, Real},  # 2D coordinates
              num = 40  :: Integer)         # Number of points for distributions 
```

```julia
foilpath = "path/to/your/airfoil.dat"   # Airfoil coordinates file path
coords = read_foil(foilpath)            # Read coordinates file
cos_foil = cosine_foil(coords, 51)      # Cosine spacing with 51 points on upper and lower surfaces
xcamthick = foil_camthick(cos_foil)     # Convert to camber-thickness representation
foiler = camthick_foil(xcamthick[:,1],  # x-components
                       xcamthick[:,2],  # Camber distribution
                       xcamthick[:,3])  # Thickness distribution
```

---

## Doublet-Source Panel Method

AeroMDAO provides convenience functions using its specific types for analyses.

The `Foil` type converts the coordinates into a friendly type for analyses. 

```julia
airfoil = naca4((2,4,1,2))
foil = Foil(airfoil)
```

The `Uniform2D` type consists of the freestream speed and angle of attack. 

```julia
V, α = 1.0, 3.0 
uniform = Uniform2D(V, α)
```

To analyse this foil with these boundary conditions using the incompressible 2D doublet-source panel method, the following method is called. Optional named arguments are provided to specify whether the source terms are non-zero, the length of the wake, and the number of panels for the analysis.

```julia
solve_case(foil                 :: Foil
           uniform              :: Uniform2D;
           sources = false      :: Bool,
           wake_length = 1e3    :: Float,
           num_panels = 60      :: Integer)
```

The method returns the lift coefficient calculated by the doublet strength of the wake panel, the lift and pressure coefficients over the panels, and the panels generated for post-processing.

<!-- ```@example aeromdao
airfoil = naca4((0,0,1,2); sharp_trailing_edge = true)  # hide
foil = Foil(airfoil)        # hide
V, α = 1.0, 5.0             # hide
uniform = Uniform2D(V, α)   # hide
cl, cls, cms, cps, panels = solve_case(foil, 
                                       uniform;
                                       viscous = false,
                                       sources = false, 
                                       wake_length = 1e3,
                                       wake_panels = 100,
                                       num_panels = 80)

println("Cl: $cl")
println("Σᵢ Clᵢ: $(sum(cls))")
println("Σᵢ Cmᵢ: $(sum(cms))")

# hide
## Pressure coefficients # hide
plot((first ∘ collocation_point).(panels), cps, # hide
     marker = 2, yflip = true, title = "Pressure Coefficient Distribution", # hide
     label = :none, xlabel = "(x/c)", ylabel = "Cp", dpi = 200) # hide
savefig("cp_plot.png"); nothing # hide
# hide
## Lift polar # hide
plot((first ∘ collocation_point).(panels), cls, # hide
     marker = 2, title = "Lift Coefficient Distribution", # hide
     label = :none, xlabel = "(x/c)", ylabel = "Cl", dpi = 200) # hide
savefig("cl_plot.png"); nothing # hide
``` -->

### Plotting

The outputs of `solve_case()` here can be used for post-processing. The following shows the visualisation setups for the pressure and lift coefficient distributions over the airfoil.

```julia
## Pressure coefficients
plot((first ∘ collocation_point).(panels), cps, 
     marker = 2, yflip = true, title = "Pressure Coefficient Distribution",
     label = :none, xlabel = "(x/c)", ylabel = "Cp")

## Lift polar
plot((first ∘ collocation_point).(panels), cls, 
     marker = 2, title = "Lift Coefficient Distribution",
     label = :none, xlabel = "(x/c)", ylabel = "Cl")
```

![cpplot](cp_plot.png)
![clplot](cl_plot.png)

---

## Wing Parametrization

The following image depicts the parametrization schema used for wing planforms in terms of nested trapezoids, consisting of airfoils and their associated chord lengths $c_n$, twist angles $\iota_n$, span lengths $b_n$, dihedrals $\delta_n$, and sweep angles $\Lambda_n$, with all angles in degrees.

![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

This parametrization is implemented as a composite type called a `HalfWing`:

```julia
HalfWing(foils      :: Vector{Foil},
         chords     :: Vector{Real},
         twists     :: Vector{Real},
         spans      :: Vector{Real},
         dihedrals  :: Vector{Real},
         sweeps     :: Vector{Real})
```

We can create a `Wing` by feeding two `HalfWing`s to it.

<!-- ```@example aeromdao
airfoil    = naca4((2,4,1,2))
foils      = Foil.(airfoil for i in 1:3)
wing_right = HalfWing(foils,  
                      [0.4, 0.2, 0.1],
                      [0., 2., 5.],
                      [1.0, 0.1], 
                      [0., 60.],  
                      [0., 30.])  

wing = Wing(wing_right, wing_right);

println("Wing —")
print_info(wing)
``` -->
---

## Vortex Lattice Method

For a 3D case, we use the vortex lattice method for initial designs, given its quick speed for fast analyses. For an analysis, you require ambient reference conditions. In this case, you need the density $\rho$ and a reference location $x_\text{ref}$ for calculating forces and moments.

Now we run the case with specifications of the number of spanwise and chordwise panels by calling the `solve_case()` function, which has an associated method:
```julia
solve_case(wing                     :: Union{Wing, HalfWing},
           freestream               :: Freestream, 
           ρ                        :: Real, 
           r_ref = [0.25, 0., 0.]   :: SVector{3, Real}; 
           span_num = 5             :: Integer, 
           chord_num = 10           :: Integer)
```

<!-- ```@example aeromdao
airfoil = naca4((2,4,1,2))  # hide
foils = Foil.(airfoil for i in 1:3) # hide
wing_right = HalfWing(foils,    # hide
                      [0.4, 0.2, 0.1],  # hide
                      [0., 2., 5.], # hide
                      [1.0, 0.1],   # hide
                      [0., 60.],    # hide
                      [0., 30.])    # hide

wing = Wing(wing_right, wing_right); # hide

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
ρ       = 1.225
ref     = [0.25 * c, 0., 0.]
U       = 10.0
α       = 2.0
β       = 2.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(U, α, β, Ω);
nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs = 
solve_case(wing, fs; 
           rho_ref = ρ, r_ref = ref, 
           area_ref = S, span_ref = b, chord_ref = c, 
           span_num = [30, 10], chord_num = 10);

print_coefficients("Wing", nf_coeffs, ff_coeffs)

nf, ff, dvs = 
solve_stability_case(wing, fs; 
                     rho_ref = ρ, r_ref = ref, 
                     area_ref = S, span_ref = b, chord_ref = c, 
                     span_num = [30, 10], chord_num = 10, name = "My Wing");

print_derivatives("Wing", dvs; farfield = true)

## Coordinates # hide
horseshoe_coords    = plot_panels(horseshoe_panels[:]) # hide
camber_coords       = plot_panels(camber_panels[:]) # hide
wing_coords         = plot_wing(wing); # hide
 # hide
CDis    = getindex.(CFs, 1) # hide
CYs     = getindex.(CFs, 2) # hide
CLs     = getindex.(CFs, 3); # hide
CL_loadings = 2sum(Γs, dims = 1)[:] / (V * b) # hide
 # hide
colpoints   = horseshoe_point.(horseshoe_panels) # hide
xs          = getindex.(colpoints, 1); # hide
ys          = getindex.(colpoints, 2); # hide
zs          = getindex.(colpoints, 3); # hide
cl_pts      = tupvector(SVector.(xs[:], ys[:], zs[:] .+ CLs[:])); # hide
 # hide
## Streamlines # hide
span_points = 20 # hide
init        = trailing_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points)  # hide
dx, dy, dz  = 0, 0, 1e-3 # hide
seed        = [ init .+ Ref([dx, dy, dz])  ; # hide
                init .+ Ref([dx, dy, -dz]) ]; # hide
 # hide
distance = 2 # hide
num_stream_points = 100 # hide
streams = plot_streams(fs, seed, horseshoes, Γs, distance, num_stream_points); # hide
 # hide
## Streamlines and panels # hide
z_limit = 1.0 # hide
plot(xaxis = "x", yaxis = "y", zaxis = "z", # hide
     aspect_ratio = 1, # hide
     camera = (45, 45), # hide
     zlim = (-0.1, z_limit), # hide
     size = (800,600)) # hide
plot!.(camber_coords,  # hide
       color = :black, label = :none) # hide
scatter!(tupvector(colpoints)[:],  # hide
         marker = 2, color = :black, label = :none) # hide
plot!.(streams,  # hide
       color = :green, label = :none) # hide
plot!(dpi = 300) # hide
savefig("streams.png"); nothing # hide
 # hide
## Span forces # hide
plot1 = plot(ys[1,:], sum(CDis, dims = 1)[:], # hide
             label = :none, xlabel = "y", ylabel = "CDi") # hide
plot2 = plot(ys[1,:], abs.(sum(CYs, dims = 1)[:]),  # hide
             label = :none, xlabel = "y", ylabel = "CY") # hide
plot3 = plot(ys[1,:], sum(CLs, dims = 1)[:],  # hide
             label = :none, xlabel = "y", ylabel = "CL") # hide
plot(plot1, plot2, plot3, layout = (3,1), size = (800,600), dpi = 300) # hide
savefig("span_forces.png"); nothing # hide
 # hide
## Lift distribution # hide
plot(xaxis = "x", yaxis = "y", zaxis = "z", # hide
     aspect_ratio = 1, # hide
     camera = (60, 45), # hide
     zlim = (-0.1, z_limit)) # hide
plot!(wing_coords, label = :none) # hide
scatter!(cl_pts, zcolor = CLs[:], label = "CL") # hide
plot!(size = (800,600), dpi = 300) # hide
 # hide
savefig("lift_dist.png"); nothing # hide
``` -->

### Plotting

The outputs of `solve_case()` here can be used for post-processing. The following shows the visualisation setups for the force distributions and streamlines over the wing.

```julia
## Coordinates
horseshoe_coords    = plot_panels(horseshoe_panels[:])
camber_coords       = plot_panels(camber_panels[:])
wing_coords         = plot_wing(wing);

CDis    = getindex.(CFs, 1)
CYs     = getindex.(CFs, 2)
CLs     = getindex.(CFs, 3);

colpoints   = horseshoe_point.(horseshoe_panels)
xs          = getindex.(colpoints, 1);
ys          = getindex.(colpoints, 2);
zs          = getindex.(colpoints, 3);
cl_pts      = tupvector(SVector.(xs[:], ys[:], zs[:] .+ CLs[:]));

## Streamlines
span_points = 20
init        = trailing_chopper(ifelse(β == 0 && Ω == zeros(3), wing.right, wing), span_points) 
dx, dy, dz  = 0, 0, 1e-3
seed        = [ init .+ Ref([dx, dy, dz])  ;
                init .+ Ref([dx, dy, -dz]) ];

distance = 2
num_stream_points = 100
streams = plot_streams(freestream, seed, horseshoes, Γs, distance, num_stream_points);

## Streamlines and panels
z_limit = 1.0
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (45, 45),
     zlim = (-0.1, z_limit),
     size = (800,600))
plot!.(camber_coords, 
       color = :black, label = :none)
scatter!(tupvector(colpoints)[:], 
         marker = 2, color = :black, label = :none)
plot!.(streams, 
       color = :green, label = :none)
plot!(dpi = 300)

## Span forces
plot1 = plot(ys[1,:], sum(CDis, dims = 1)[:],
             label = :none, xlabel = "y", ylabel = "CDi")
plot2 = plot(ys[1,:], abs.(sum(CYs, dims = 1)[:]), 
             label = :none, xlabel = "y", ylabel = "CY")
plot3 = plot(ys[1,:], sum(CLs, dims = 1)[:], 
             label = :none, xlabel = "y", ylabel = "CL")
plot(plot1, plot2, plot3, layout = (3,1), size = (800,600), dpi = 300)

## Lift distribution
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     aspect_ratio = 1,
     camera = (60, 45),
     zlim = (-0.1, z_limit))
plot!(wing_coords, label = :none)
scatter!(cl_pts, zcolor = CLs[:], label = "CL")
plot!(size = (800,600), dpi = 300)
```

![streams](streams.png)
![spanforces](span_forces.png)
![liftdist](lift_dist.png)