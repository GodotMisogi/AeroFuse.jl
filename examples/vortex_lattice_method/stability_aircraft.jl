## Aircraft stability analysis example
using AeroMDAO
using LinearAlgebra # For norm()

## Lifting surfaces

# Wing
wing  = WingSection(span       = 8.0,
                    dihedral   = 5.0,
                    sweep_LE   = 15.0,
                    taper      = 0.4,
                    root_chord = 2.0,
                    root_twist = 0.0,
                    tip_twist  = -2.0,
                    root_foil  = naca4((2,4,1,2)),
                    tip_foil   = naca4((2,4,1,2)),
                    position   = [0., 0., 0.])

# Horiontal tail
htail = WingSection(span       = 2.0,
                    dihedral   = 0.0,
                    sweep_LE   = 15.0,
                    taper      = 0.6,
                    root_chord = 0.8,
                    root_twist = 0.0,
                    tip_twist  = 0.0,
                    root_foil  = naca4((0,0,1,2)),
                    tip_foil   = naca4((0,0,0,9)),
                    position   = [5., 0., 0.],
                    angle      = 0.,
                    axis       = [0., 1., 0.]);

# Vertical tail
vtail = HalfWingSection(span       = 0.8,
                        dihedral   = 0.0,
                        sweep_LE   = 8.0,
                        taper      = 0.6,
                        root_chord = 0.8,
                        root_twist = 0.0,
                        tip_twist  = 0.,
                        root_foil  = naca4((0,0,0,9)),
                        tip_foil   = naca4((0,0,0,9)),
                        position   = [5., 0., 0.],
                        angle      = 90.,
                        axis       = [1., 0., 0.])
vtail_plan = plot_wing(vtail)           

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## Static stability
wing_mac  = mean_aerodynamic_center(wing)
htail_mac = mean_aerodynamic_center(htail)
vtail_mac = mean_aerodynamic_center(vtail)

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
x_w = [ wing_mac[1], 0, 0 ]

# Horizontal tail volume coefficient
S_h = projected_area(htail)
x_h = [ htail_mac[1], 0, 0 ]
l_h = norm(x_h - x_w)
V_h = S_h * l_h / (S * c)
println("Horizontal TVC      V_h: $V_h")

# Vertical tail volume coefficient
S_v = projected_area(vtail)
x_v = [ vtail_mac[1], 0, 0 ]
l_v = norm(x_v - x_w)
V_v = S_v * l_v / (S * b)
println("Vertical TVC        V_v: $V_v")

## Panelling and assembly
wing_panels, wing_normals   = panel_wing(wing, [20], 10)
htail_panels, htail_normals = panel_wing(htail, [10], 5)
vtail_panels, vtail_normals = panel_wing(vtail, [10], 5)

aircraft = Dict("Wing"            => Horseshoe.(wing_panels,  wing_normals ),
                "Horizontal Tail" => Horseshoe.(htail_panels, htail_normals),
                "Vertical Tail"   => Horseshoe.(vtail_panels, vtail_normals))

## Evaluate case
ac_name = "My Aircraft"
ρ 		= 1.225
ref     = x_w
V, α, β = 1.0, 3.0, 0.0
Ω 		= [0.0, 0.0, 0.0]
fs 	    = Freestream(V, α, β, Ω)

@time dv_data = 
solve_stability_case(aircraft, fs; 
                     rho_ref     = ρ,           # Reference density
                     r_ref       = ref,         # Reference point for moments
                     area_ref    = S,           # Reference area
                     span_ref    = b,           # Reference span
                     chord_ref   = c,           # Reference chord
                     name        = ac_name,     # Aircraft name
                     print       = true,        # Prints the results for the entire aircraft
                     print_components = true,   # Prints the results for each component
                    );

## Process data
labels = (collect ∘ keys)(dv_data) # Get aircraft component names from analysis
comp = labels[1]                   # Pick your component
nf, ff, dvs = dv_data[comp];       # Get the nearfield, farfield, and stablity derivative coefficients
print_case(dv_data, comp)          # Pretty-print the results

## Aerodynamic quantities of aircraft
nf_plane, ff_plane, dvs_plane = dv_data[labels[1]]

# Center of pressure
Cm   = nf_plane[5]	# Moment coefficient
CL   = nf_plane[3]	# Lift coefficient
cp   = -c * Cm / CL

# Neutral point
Cm_α = dvs_plane[5,1]
CL_α = dvs_plane[3,1]
np   = -c * Cm_α / CL_α

# Spiral stability
Cl_β = dvs_plane[4,2]
Cl_r = dvs_plane[4,5]
Cn_r = dvs_plane[6,5]
Cn_β = dvs_plane[6,2]

γ = Cl_β * Cn_r / (Cl_r * Cn_β) # Check degree-radian issue

# Locations
x_np = [ ref[1] .+ np; zeros(2) ]	# Neutral point
x_cp = [ ref[1] .+ cp; zeros(2) ]	# Center of pressure

println("Aerodynamic Center x_ac: $(x_w[1]) m")
println("Neutral Point      x_np: $(x_np[1]) m")
println("Center of Pressure x_cp: $(x_cp[1]) m")
println("Spiral Stability      γ: $γ")

## Plotting everything
using Plots
pyplot()

wing_plan = plot_wing(wing)

#
z_limit = b
aircraft_plot = 
plot(xaxis = "x", yaxis = "y", zaxis = "z",
     camera = (30, 60),
     xlim = (0, z_limit),
    #  ylim = (-z_limit/2, z_limit/2),
     zlim = (-z_limit/2, z_limit/2),
     size = (1920, 1080),
    #  legend = false
    )

plot!(wing_plan, label = "Wing")
plot!(htail_plan, label = "Horizontal Tail")
plot!(vtail_plan, label = "Vertical Tail")

scatter!(Tuple(wing_mac), color = :blue, label = "Wing MAC")
scatter!(Tuple(htail_mac), color = :red, label = "Horizontal Tail MAC")
scatter!(Tuple(vtail_mac), color = :green, label = "Vertical Tail MAC")

scatter!(Tuple(x_np), color = :orange, label = "Neutral Point")
scatter!(Tuple(x_cp), color = :brown, label = "Center of Pressure")
# savefig(aircraft_plot, "Aircraft.png")
plot!()