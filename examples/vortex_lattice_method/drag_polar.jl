## Drag polar analysis case
using AeroMDAO
using DataFrames
using Plots

## Lifting surfaces

# Wing
wing_foils = fill(naca4((0,0,1,2)), 2)
wing       = Wing(foils     = wing_foils,
                  chords    = [1.0, 0.6],
                  twists    = [0.0, 0.0],
                  spans     = [5.0],
                  dihedrals = [0.],
                  LE_sweeps = [2.29]);

# Horizontal tail
htail_foils = fill(naca4((0,0,1,2)), 2)
htail       = Wing(foils     = htail_foils,
                   chords    = [0.7, 0.42],
                   twists    = [0.0, 0.0],
                   spans     = [1.25],
                   dihedrals = [0.],
                   LE_sweeps = [6.39],
                   position	 = [4., 0, 0],
                   angle     = deg2rad(-0.),
                   axis      = [0., 1., 0.])


# Vertical tail
vtail_foils = fill(naca4((0,0,1,2)), 2)
vtail = HalfWing(foils     = vtail_foils,
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 LE_sweeps = [7.97],
                 position  = [4., 0, 0],
                 angle     = 90.,
                 axis      = [1., 0., 0.])

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

# Assembly
wing_panels, wing_normals   = panel_wing(wing, [20], 10)
htail_panels, htail_normals = panel_wing(htail, [12], 12)
vtail_panels, vtail_normals = panel_wing(vtail, [12], 10)

aircraft = ComponentArray(wing  = Horseshoe.(wing_panels,  wing_normals ),
                          htail = Horseshoe.(htail_panels, htail_normals),
                          vtail = Horseshoe.(vtail_panels, vtail_normals))


S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

function vlm_analysis(aircraft, fs, refs, print = false)
    # Evaluate case
    system = solve_case(aircraft, fs, refs;
                        print = print);

    # Get data
    [ farfield(system); nearfield(system) ]
end

##
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, β    = 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
αs      = -5:0.5:5
fses    = map(α -> Freestream(V, α, β, Ω), αs);
refs    = References(ρ, S, b, c, ref)

## Evaluate cases
results = vlm_analysis.(Ref(aircraft), fses, Ref(refs), true);

## Generate DataFrame
data = DataFrame((α, res...) for (α, res) in zip(αs, results))
rename!(data, [:α, :CD_ff, :CY_ff, :CL_ff, :CD_nf, :CY_nf, :CL_nf, :Cl, :Cm, :Cn])

##
# pyplot(dpi = 300)
plot(xlabel = "CD", ylabel = "CL", title = "Drag Polar")
plot!(data[!,"CD_ff"], data[!,"CL_ff"], label = :none, marker = :dot)