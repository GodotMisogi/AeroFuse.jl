## Example script for building surrogate models using VLM analysis, needs Surrogates and DataFrames installed
using Surrogates
using DataFrames
using AeroMDAO

## Lifting surfaces

# Wing
wing_foils = Foil.(fill(naca4((0,0,1,2)), 2))
wing       = Wing(foils     = wing_foils,
                  chords    = [1.0, 0.6],
                  twists    = [0.0, 0.0],
                  spans     = [5.0],
                  dihedrals = [0.],
                  sweep_LEs = [2.29]);

# Horizontal tail
htail_foils = Foil.(fill(naca4((0,0,1,2)), 2))
htail       = Wing(foils     = htail_foils,
                   chords    = [0.7, 0.42],
                   twists    = [0.0, 0.0],
                   spans     = [1.25],
                   dihedrals = [0.],
                   sweep_LEs = [6.39],
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
                 sweep_LEs = [7.97],
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
wing_mac = mean_aerodynamic_center(wing)

## VLM setup
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

## Evaluate one case for compilation and test

# Case
ac_name = "My Aircraft"
ρ       = 1.225
ref     = [ wing_mac[1], 0, 0 ]
V, α, β = 1.0, 0.0, 0.0
Ω       = [0.0, 0.0, 0.0]
fs      = Freestream(V, α, β, Ω)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

@time coeffs = vlm_analysis(aircraft, fs, ρ, ref, S, b, c, true)

## Looping
αs, βs   = float.(-4:4), float.(-3:3)
fses     = [ Freestream(V, α, β, Ω) for α in αs, β in βs ]

## Evaluate cases
results = vlm_analysis.(Ref(aircraft), fses, ρ, Ref(ref), S, b, c, true);

## Generate DataFrame
αβs  = [ (α, β) for α in αs, β in βs ]
data = DataFrame([ (α, β, CD, CY, CL, Cl, Cm, Cn) for ((α, β), (CD, CY, CL, Cl, Cm, Cn)) in zip(αβs, results) ][:])
rename!(data, [:α, :β, :CD, :CY, :CL, :Cl, :Cm, :Cn])

##
using Plots
pyplot(size = (1280, 720), markersize = 2, dpi = 300)

## The labels don't work for some reason, need to fix that later...
CD_plot = scatter(data[:,1], data[:,2], data[:,3], label = "CD")
CY_plot = scatter(data[:,1], data[:,2], data[:,4], label = "CY")
CL_plot = scatter(data[:,1], data[:,2], data[:,5], label = "CL")
Cl_plot = scatter(data[:,1], data[:,2], data[:,6], label = "Cl")
Cm_plot = scatter(data[:,1], data[:,2], data[:,7], label = "Cm")
Cn_plot = scatter(data[:,1], data[:,2], data[:,8], label = "Cn")

plots = [ CD_plot CY_plot CL_plot Cl_plot Cm_plot Cn_plot ]
plot(plots..., layout = (2, 3), xlabel = "α, ᵒ", ylabel = "β, ᵒ")
plot!()

## Surrogate training data
lower_bound, upper_bound = minimum.(eachcol(data[:,3])), maximum.(eachcol(data[:,3]))
zs = data[:,3]

## Generate surrogate
surrogate_model = LobachevskySurrogate(αβs, zs, lower_bound, upper_bound);

## Evaluate surrogate
num     = 100
xs, ys  = range(-4, 4, length = num), range(-3, 3, length = num)
xys     = [ (x, y) for x in xs, y in ys ]
surr_zs = surrogate_model.(xys)

## Create DataFrame
surr_data = DataFrame((x, y, z) for ((x, y), z) in zip(xys, surr_zs))
rename!(surr_data, [:α, :β, :CD])

## Plot surrogate
scatter(data[:,1], data[:,2], data[:,3], xlabel = "α, ᵒ", ylabel = "β, ᵒ", zlabel = "CD", label = "CD")
mesh3d!(surr_data[:,1], surr_data[:,2], surr_data[:,3], label = "Surrogate CD")