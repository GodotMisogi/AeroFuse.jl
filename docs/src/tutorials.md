```@meta
EditURL = "<unknown>/tutorials.jl"
```

```@setup aeromdao
using AeroMDAO
using Plots
using LaTeXStrings
plotlyjs()
```


# Tutorials

## Your First Airfoil

### Plotting

### Aerodynamic Analysis

## Your First Wing

### Parametrization

!!! definition
    A **wing section** consists of two foil profiles and their chord lengths and twist angles. Between them is their span length with associated _leading-edge_ dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for spans, dihedrals, sweeps, and ``n+1`` entries for foils, chords, and twists for some ``n \in \mathbb N``
![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)
You can define a wing in this parametrization by using the `HalfWing` type:
```julia
 HalfWing(foils          :: Vector{Foil}, # Foil profiles
          chords         :: Vector{Real}, # Chord lengths
          twists         :: Vector{Real}, # Twist angles
          spans          :: Vector{Real}, # Section span lengths
          dihedrals      :: Vector{Real}, # Dihedral angles
          LE_sweeps      :: Vector{Real}) # Leading-edge sweep angles
```

First, we define a `Vector` of `Foil`s.

```@example tutorials
airfoil_1 = Foil(naca4(4,4,1,2))
airfoils  = [ airfoil_1, Foil(naca4(0,0,1,2)) ]
```

Now that we have our foil profiles in the appropriate type, we can define our wing:"

```@example tutorials
wing_right = Wing(foils     = airfoils,
                  chords    = [1.0, 0.6],
                  twists    = [2.0, 0.0],
                  spans     = [4.0],
                  dihedrals = [5.],
                  LE_sweeps = [5.])

chords(wing)
```

### Visualization

Now let's see what the outline of our wing looks like, using the following function to get the coordinates.

```@example tutorials
wing_outline = plot_wing(wing)

plt1 = plot(
            xaxis = "x", yaxis = "y", zaxis = "z",
            aspect_ratio = 1,
            camera = (30, 45),
            zlim = (-0.1, span(wing) / 2),
           )
plot!(wing_outline[:,1], wing_outline[:,2], wing_outline[:,3], label = "Wing")
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

