## Aircraft stability analysis example
using AeroFuse

## Lifting surfaces

# Wing
wing  = WingSection(
    area       = 6.4,
    aspect     = 10.,
    dihedral   = 5.0,
    sweep      = 15.0,
    taper      = 0.4,
    root_twist = 0.0,
    tip_twist  = -2.0,
    root_foil  = naca4(2,4,1,2),
    tip_foil   = naca4(2,4,0,9),
    position   = [0., 0., 0.],
    symmetry   = true
)

# Horizontal tail
htail = WingSection(
    area       = 2.56 / 2,
    aspect     = 6.25,
    dihedral   = 0.0,
    sweep      = 15.0,
    taper      = 0.6,
    root_twist = 0.0,
    tip_twist  = 0.0,
    root_foil  = naca4(0,0,1,2),
    tip_foil   = naca4(0,0,0,9),
    position   = [5., 0., -0.1],
    angle      = 0.,
    axis       = [0., 1., 0.],
    symmetry   = true
);

# Vertical tail
vtail = WingSection(
    area       = 0.512 / 2,
    aspect     = 1.25,
    dihedral   = 0.0,
    sweep      = 8.0,
    taper      = 0.6,
    root_twist = 0.0,
    tip_twist  = 0.0,
    root_foil  = naca4(0,0,0,9),
    tip_foil   = naca4(0,0,0,9),
    position   = [5., 0., 0.],
    angle      = 90.,
    axis       = [1., 0., 0.]
)

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
fs = Freestream(
    alpha = 1.0, 
    beta  = 0.0, 
    omega = [0.,0.,0.]
)
refs = References(
    speed    = 10.0,
    density  = 1.225, 
    area     = projected_area(wing),
    span     = span(wing), 
    chord    = mean_aerodynamic_chord(wing),
    location = mean_aerodynamic_center(wing)
)

ac_name = :aircraft
@time dv_data = freestream_derivatives(
    aircraft, fs, refs;
#    axes             = Wind(),
    name             = ac_name,
    print            = true,    # Prints the results for only the aircraft
    # print_components = true,    # Prints the results for all components
);

## Aerodynamic quantities of aircraft
ac_dvs = dv_data[ac_name]

# Longitudinal stability

## Center of pressure
x_cp = -refs.chord * ac_dvs.Cm / ac_dvs.CZ

## Neutral point
x_np = -refs.chord * ac_dvs.Cm_al / ac_dvs.CZ_al

## Spiral stability
γ = ac_dvs.Cl_be * ac_dvs.Cn_rb / (ac_dvs.Cl_rb * ac_dvs.Cn_be)

## Other quantities
htail_dvs = dv_data.htail

# Locations
x_np, y_np, z_np = loc_np = [ refs.location[1] .+ x_np; zeros(2) ]  # Neutral point
x_cp, y_cp, z_cp = loc_cp = [ refs.location[1] .+ x_cp; zeros(2) ]  # Center of pressure

println("Aerodynamic Center x_ac: $(x_w) m")
println("Neutral Point      x_np: $(x_np) m")
println("Center of Pressure x_cp: $(x_cp) m")
println("Spiral Stability      γ: $γ")

## Plotting everything
using Plots

gr(dpi = 300)

#
z_limit = span(wing)
aircraft_plot = plot(
    xaxis = "x", yaxis = "y", zaxis = "z",
    camera = (30, 30),
    xlim = (0, z_limit),
    # ylim = (-z_limit/2, z_limit/2),
    zlim = (-z_limit/2, z_limit/2),
    size = (900, 600),
    # legend = false
    )

plot!(wing,  label = "Wing")
plot!(htail, label = "Horizontal Tail")
plot!(vtail, label = "Vertical Tail")

scatter!(Tuple(loc_np), color = :orange, label = "Neutral Point")
scatter!(Tuple(loc_cp), color = :brown, label = "Center of Pressure")
savefig(aircraft_plot, "plots/static_stability.png")
plot!()


## Alpha sweep
function alpha_sweep(aircraft, refs, α)
    ## Evaluate case
    fs      = Freestream(alpha = α, 
                         beta  = 0.0, 
                         omega = [0.,0.,0.])

    system = solve_case(aircraft, fs, refs)

    dvs = freestream_derivatives(system)

    return dvs
end

##
αs  = -8:0.5:8
@time res = map(α -> alpha_sweep(aircraft, refs, α), αs);

## Data processing

## Lift coefficient
CLs = [ re.aircraft.CZ for re in res ]

plot(αs, CLs, xlabel = "α", ylabel = "CL")

## Moment coefficient
Cms = [ re.aircraft.Cm for re in res ]

plot(αs, Cms, xlabel = "α", ylabel = "Cm")

# Cm-CL (only relevant for non-linear CL-α variation)
# plot(CLs, Cms, xlabel = "CL", ylabel = "Cm")

##
Cm_αs  = [ re.aircraft.Cm_al for re in res ]
CL_αs  = [ re.aircraft.CZ_al for re in res ]

##
plot(CL_αs, Cm_αs, xlabel = "CL_α", ylabel = "Cm_α")

##
∂Cm_∂CLs = Cm_αs ./ CL_αs

plot(αs, ∂Cm_∂CLs, ls = :solid, xlabel = "α", ylabel = "Cm/CL")