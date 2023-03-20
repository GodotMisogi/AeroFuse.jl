## Wing analysis case
using AeroFuse

## Wing section setup
wing = Wing(
    foils     = fill(naca4((2,4,1,2)), 3),
    chords    = [2.0, 1.2, 0.6],
    twists    = [0.0, 0.0, 0.0],
    spans     = [5., 0.6],
    dihedrals = [5., 5.],
    sweeps    = [20.,30.],
    w_sweep   = 0.25, # Quarter-chord sweep
    symmetry  = true,
    position  = [4.0, 0.0, 0.0]
    # flip      = true
)

##
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Meshing and assembly
wing_mesh = WingMesh(wing, [24,6], 6, 
                    #  span_spacing = Cosine()
                    );
aircraft = ComponentVector(wing = make_horseshoes(wing_mesh))

# Freestream conditions
fs  = Freestream(
    alpha = 2.0, # deg
    beta  = 0.0, # deg
    omega = [0.,0.,0.]
)

# Reference values
ref = References(
    speed     = 150, # m/s
    density   = 1.225, # kg/m³
    viscosity = 1.5e-5, # ???
    area      = projected_area(wing), # m²
    span      = span(wing), # m
    chord     = mean_aerodynamic_chord(wing), # m
    location  = mean_aerodynamic_center(wing) # m
)

## Solve system
@time begin 
    system = solve_case(
        aircraft, fs, ref;
        compressible = true,
        # print            = true, # Prints the results for only the aircraft
        # print_components = true, # Prints the results for all components
    );

    # Compute dynamics
    ax       = Wind() # Geometry, Stability(), Body()
    CFs, CMs = surface_coefficients(system; axes = ax)
    # Fs, Ms   = surface_dynamics(system; axes = ax)
    # Fs       = surface_forces(system; axes = ax)
    # vels     = surface_velocities(system)

    nf  = nearfield(system)
    ff  = farfield(system)

    nfs = nearfield_coefficients(system)
    ffs = farfield_coefficients(system)
        
    # Force/moment coefficients and derivatives
    dvs = freestream_derivatives(system; 
        axes = ax,
        print_components = true,
        farfield = true
    )
end;

## Viscous drag prediction using empirical models

# Equivalent flat-plate skin-friction estimation
x_tr        = fill(0.98, 2) # Transition locations over sections
CDv_plate   = profile_drag_coefficient(wing_mesh, x_tr, ref)

## Local-dissipation drag estimation
cam_panels  = camber_panels(wing_mesh)
edge_speeds = surface_velocities(system).wing
CDv_diss    = profile_drag_coefficient(wing_mesh, x_tr, edge_speeds, ref)

## Viscous drag coefficient
CDv = CDv_diss

## Total force coefficients with viscous drag prediction
CDi_nf, CY_nf, CL_nf, Cl, Cm, Cn = nf = nearfield(system) 
CDi_ff, CY_ff, CL_ff = ff = farfield(system)

nf_v = [ CDi_nf + CDv; CDv; nf ]
ff_v = [ CDi_ff + CDv; CDv; ff ]

print_coefficients(nf_v, ff_v)

## Plotting
using Plots
gr()

## Coordinates
Plots.plot(
    aspect_ratio = 1,
    camera = (30, 30),
    zlim = span(wing) .* (-0.5, 0.5),
    size = (800, 600)
)
Plots.plot!(wing_mesh, label = "Wing")
Plots.plot!(system, wing, dist = 3, num_stream = 50, span = 10, color = :green)

## Compute spanwise loads
span_loads = spanwise_loading(wing_mesh, ref, CFs.wing, system.circulations.wing)

## Plot spanwise loadings
plot_CD = plot(span_loads[:,1], span_loads[:,2], label = :none, ylabel = "CDi")
plot_CY = plot(span_loads[:,1], span_loads[:,3], label = :none, ylabel = "CY")
plot_CL = begin
            plot(span_loads[:,1], span_loads[:,4], label = :none, xlabel = "y", ylabel = "CL")
            plot!(span_loads[:,1], span_loads[:,5], label = "Normalized", xlabel = "y")
          end
plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))

## VARIABLE ANALYSES
#=========================================================#

using Accessors
using Base.Iterators: product

## Speed sweep
Vs = 1.0:10:300
res_Vs = combinedimsview(
    map(Vs) do V
        refs = @set ref.speed = V
        sys = solve_case(aircraft, fs, refs)
        [ mach_number(refs); farfield(sys)...; nearfield(sys)... ]
    end, (1))

plot(
    res_Vs[:,1], res_Vs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "M",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)

## Alpha sweep
αs = -5:0.5:5
res_αs = combinedimsview(
    map(αs) do α
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref)
        [ α; farfield(sys); nearfield(sys) ]
    end, (1)
)

plot(
    res_αs[:,1], res_αs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "α",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)

## (Speed, alpha) sweep
res = combinedimsview(
    map(product(Vs, αs)) do (V, α)
        refs = @set ref.speed = V
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, refs)
        [ mach_number(refs); α; farfield(sys); nearfield(sys) ]
    end
)

res_p = permutedims(res, (3,1,2))

## CDi
plt_CDi_ff = plot(camera = (60,45))
[ plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,3,n], 
    ylabel = "α", xlabel = "M", zlabel = "CDi_ff", 
    label = "", c = :black,
) for n in axes(res_p,3) ]

# CL
plt_CL_ff = plot(camera = (45,45))
[ plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,8,n], 
    ylabel = "α", xlabel = "M", zlabel = "CL_ff", 
    label = "", c = :black,
) for n in axes(res_p,3) ]

plt_Cm_ff = plot(camera = (45,45))
[ plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,10,n], 
    ylabel = "α", xlabel = "M", zlabel = "Cm", 
    label = "", c = :black,
) for n in axes(res_p,3) ]

##
plot(plt_CDi_ff, plt_CL_ff, plt_Cm_ff, layout = (1,3), size = (1300, 400))
