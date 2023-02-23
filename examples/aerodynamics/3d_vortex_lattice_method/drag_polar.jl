## Drag polar analysis case
using AeroFuse
using DataFrames
using Plots

## Lifting surfaces

## Surfaces

# Wing
wing = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [1.0, 0.6],
    twists    = [0.0, 0.0],
    spans     = [4.0],
    dihedrals = [0.],
    sweeps    = [5.],
    symmetry  = true,
);

# Horizontal tail
htail = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25],
    dihedrals = [0.],
    sweeps    = [6.39],
    position  = [4., 0, 0.],
    angle     = 0.,
    axis      = [0., 1., 0.],
    symmetry  = true
)

# Vertical tail
vtail = Wing(
    foils     = fill(naca4(0,0,0,9), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [0.5],
    dihedrals = [0.],
    sweeps    = [7.97],
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.],
)

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## WingMesh type
wing_mesh  = WingMesh(wing, [12], 6, 
                span_spacing = Cosine()
            )
htail_mesh = WingMesh(htail, [12], 6, 
                span_spacing = Cosine()
            )
vtail_mesh = WingMesh(vtail, [12], 6, 
                span_spacing = Cosine()
            )

aircraft = ComponentArray(
    wing  = make_horseshoes(wing_mesh),
    htail = make_horseshoes(htail_mesh),
    vtail = make_horseshoes(vtail_mesh)
);

## Case
fs  = Freestream(
    alpha = 0.0, 
    beta  = 0.0, 
    omega = [0., 0., 0.]
);

ref = References(
    speed     = 1.0, 
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing)
)

## Evaluate cases
using Accessors

αs = -5:0.5:5
results = combinedimsview(
    map(αs) do α
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref)
        [ α; farfield(sys)...; nearfield(sys)... ]
    end, (1)
)

## Generate DataFrame
data = DataFrame(results, [:α, :CD_ff, :CY_ff, :CL_ff, :CD_nf, :CY_nf, :CL_nf, :Cl, :Cm, :Cn])

##
# pyplot(dpi = 300)
plot(xlabel = "CD", ylabel = "CL", title = "Drag Polar")
plot!(data[!,"CD_ff"], data[!,"CL_ff"], label = :none, ms = :dot)