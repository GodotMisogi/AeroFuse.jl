## Drag polar analysis case
using AeroMDAO
using DataFrames
using Plots

## Lifting surfaces

# Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing       = Wing(foils     = wing_foils,
                  chords    = [1.0, 0.6],
                  twists    = [0.0, 0.0],
                  spans     = [5.0],
                  dihedrals = [0.],
                  LE_sweeps = [2.29]);

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
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
vtail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
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

aircraft = Dict("Wing"            => Horseshoe.(wing_panels,  wing_normals ),
                "Horizontal Tail" => Horseshoe.(htail_panels, htail_normals),
                "Vertical Tail"   => Horseshoe.(vtail_panels, vtail_normals))


S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
    # Evaluate case
    data =  solve_case(aircraft, fs;
                       rho_ref   = ρ,
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print     = print
                      );

    # Get data
    nf_coeffs, ff_coeffs = data["Aircraft"][1:2]

    # Filter relevant data
    [ ff_coeffs; nf_coeffs[4:6] ]
end

##
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [0.25c, 0., 0.]
V, β    = 1.0, 0.0
Ω       = [0.0, 0.0, 0.0]
αs      = -5:0.5:5
fses    = map(α -> Freestream(V, α, β, Ω), αs);

## Evaluate cases
results = vlm_analysis.(Ref(aircraft), fses, ρ, Ref(ref), S, b, c, true);

## Generate DataFrame
data = DataFrame([ (α, CD, CY, CL, Cl, Cm, Cn) for (α, (CD, CY, CL, Cl, Cm, Cn)) in zip(αs, results) ][:])
rename!(data, [:α, :CD, :CY, :CL, :Cl, :Cm, :Cn])

##
pyplot(dpi = 300)
plot(xlabel = "CD", ylabel = "CL", title = "Drag Polar")
plot!(data[!,"CD"], data[!,"CL"], label = :none, marker = :dot)