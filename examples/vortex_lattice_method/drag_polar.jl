## Drag polar analysis case
using AeroMDAO
using DataFrames
using Plots

## Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing       = Wing(foils     = wing_foils,
                  chords    = [1.0, 0.6],
                  twists    = [0.0, 0.0],
                  spans     = [5.0],
                  dihedrals = [11.3],
                  sweep_LEs = [2.29]);
print_info(wing, "Wing")

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
htail       = Wing(foils     = htail_foils,
                   chords    = [0.7, 0.42],
                   twists    = [0.0, 0.0],
                   spans     = [1.25],
                   dihedrals = [0.],
                   sweep_LEs = [6.39])
print_info(htail, "Horizontal Tail")

# Vertical tail
vtail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
vtail = HalfWing(foils     = vtail_foils, 
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.0],
                 dihedrals = [0.],
                 sweep_LEs = [7.97])
print_info(vtail, "Vertical Tail")

# Assembly
wing_panels  = panel_wing(wing, [20], 10);
htail_panels = panel_wing(htail, [12], 12;
                          position	= [4., 0, 0],
                          angle 	= deg2rad(-0.),
                          axis 	  	= [0., 1., 0.]
                         )
vtail_panels = panel_wing(vtail, [12], 10; 
                          position 	= [4., 0, 0],
                          angle 	= π/2, 
                          axis 	 	= [1., 0., 0.]
                         )

aircraft = Dict(
                "Wing" 			  	=> wing_panels,
                # "Horizontal Tail" 	=> htail_panels,
                # "Vertical Tail"   	=> vtail_panels
                )

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
    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horses, Γs = data["Aircraft"]

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
plotly(dpi = 300)
plot(xlabel = "CD", ylabel = "CL", title = "Drag Polar")
plot!(data[!,"CD"], data[!,"CL"], label = :none, marker = :dot)