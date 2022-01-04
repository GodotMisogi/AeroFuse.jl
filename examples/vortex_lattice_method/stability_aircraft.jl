## Aircraft stability analysis example
using AeroMDAO

## Lifting surfaces

# Wing
wing  = WingSection(span       = 8.0,
                    dihedral   = 5.0,
                    sweep_LE   = 15.0,
                    taper      = 0.4,
                    root_chord = 2.0,
                    root_twist = 0.0,
                    tip_twist  = -2.0,
                    root_foil  = naca4(2,4,1,2),
                    tip_foil   = naca4(2,4,1,2),
                    position   = [0., 0., 0.])

# Horiontal tail
htail = WingSection(span       = 2.0,
                    dihedral   = 0.0,
                    sweep_LE   = 15.0,
                    taper      = 0.6,
                    root_chord = 0.8,
                    root_twist = 0.0,
                    tip_twist  = 0.0,
                    root_foil  = naca4(0,0,1,2),
                    tip_foil   = naca4(0,0,0,9),
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
                        root_foil  = naca4(0,0,0,9),
                        tip_foil   = naca4(0,0,0,9),
                        position   = [5., 0., 0.],
                        angle      = 90.,
                        axis       = [1., 0., 0.])

# Print info
print_info(wing, "Wing")
print_info(htail, "Horizontal Tail")
print_info(vtail, "Vertical Tail")

## Static stability
x_w, y_w, z_w = wing_mac  = mean_aerodynamic_center(wing)
x_h, y_h, z_h = htail_mac = mean_aerodynamic_center(htail)
x_v, y_v, z_v = vtail_mac = mean_aerodynamic_center(vtail)

S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)

# Horizontal tail volume coefficient
S_h = projected_area(htail)
l_h = x_h - x_w
V_h = S_h * l_h / (S * c)
println("Horizontal TVC      V_h: $V_h")

# Vertical tail volume coefficient
S_v = projected_area(vtail)
l_v = x_v - x_w
V_v = S_v * l_v / (S * b)
println("Vertical TVC        V_v: $V_v")

## Meshing and assembly
wing_mesh  = WingMesh(wing, [20], 10)
htail_mesh = WingMesh(htail, [10], 5)
vtail_mesh = WingMesh(vtail, [10], 5)

aircraft = ComponentArray(
                          wing  = make_horseshoes(wing_mesh),
                          htail = make_horseshoes(htail_mesh),
                          vtail = make_horseshoes(vtail_mesh)
                         );

## Evaluate case
ac_name = :aircraft
fs      = Freestream(alpha = 1.0, 
                     beta  = 0.0, 
                     omega = [0.,0.,0.])

# Reference values
refs    = References(
                     speed    = 10.0,
                     density  = 1.225, 
                     area     = projected_area(wing),   
                     span     = span(wing), 
                     chord    = mean_aerodynamic_chord(wing), 
                     location = mean_aerodynamic_center(wing)
                    )

@time dv_data =
solve_case_derivatives(aircraft, fs, refs;
                     print            = true,    # Prints the results for only the aircraft
                     print_components = true,    # Prints the results for all components
                    );

## Aerodynamic quantities of aircraft
nf_plane  = dv_data.aircraft.NF
ff_plane  = dv_data.aircraft.NF
dvs_plane = dv_data.aircraft.dNF

# Center of pressure
Cm   = nf_plane[5]  # Moment coefficient
CL   = nf_plane[3]  # Lift coefficient
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

γ    = Cl_β * Cn_r / (Cl_r * Cn_β) # Check degree-radian issue

## Other quantities
dvs_htail = dv_data.htail.dNF
CL_α_h    = dvs_htail[3,1]

# Locations
x_np = [ ref[1] .+ np; zeros(2) ]  # Neutral point
x_cp = [ ref[1] .+ cp; zeros(2) ]  # Center of pressure

println("Aerodynamic Center x_ac: $(x_w[1]) m")
println("Neutral Point      x_np: $(x_np[1]) m")
println("Center of Pressure x_cp: $(x_cp[1]) m")
println("Spiral Stability      γ: $γ")

## Plotting everything
using Plots

wing_plan  = plot_wing(wing)
htail_plan = plot_wing(htail)
vtail_plan = plot_wing(vtail)

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

plot!(wing_plan,  label = "Wing")
plot!(htail_plan, label = "Horizontal Tail")
plot!(vtail_plan, label = "Vertical Tail")

scatter!(Tuple(wing_mac),  color = :blue,  label = "Wing MAC")
scatter!(Tuple(htail_mac), color = :red,   label = "Horizontal Tail MAC")
scatter!(Tuple(vtail_mac), color = :green, label = "Vertical Tail MAC")

scatter!(Tuple(x_np), color = :orange, label = "Neutral Point")
scatter!(Tuple(x_cp), color = :brown, label = "Center of Pressure")
# savefig(aircraft_plot, "Aircraft.png")
plot!()