## Wing analysis case
using AeroFuse

## Surfaces

# Wing
wing = Wing(
    foils     = fill(naca4(2,4,1,2), 2),
    chords    = [1.0, 0.6],
    twists    = [2.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.],
    sweeps    = [5.],
    symmetry  = true,
    # flip      = true
);

x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing);

# Horizontal tail
htail = Wing(
    foils     = fill(naca4(0,0,1,2), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.25],
    dihedrals = [0.],
    sweeps    = [6.39],
    position  = [4., 0, -0.1],
    angle     = -2.,
    axis      = [0., 1., 0.],
    symmetry  = true
)

# Vertical tail
vtail = Wing(
    foils     = fill(naca4(0,0,0,9), 2),
    chords    = [0.7, 0.42],
    twists    = [0.0, 0.0],
    spans     = [1.0],
    dihedrals = [0.],
    sweeps    = [7.97],
    position  = [4., 0, 0],
    angle     = 90.,
    axis      = [1., 0., 0.]
)

## WingMesh type
wing_mesh = WingMesh(wing, [24], 6,
    # span_spacing = Cosine()
)
htail_mesh = WingMesh(htail, [12], 6, 
    # span_spacing = Cosine()
)
vtail_mesh = WingMesh(vtail, [12], 6, 
    # span_spacing = Cosine()
)

aircraft =  ComponentVector(
    wing  = make_vortex_rings(wing_mesh),
    htail = make_vortex_rings(htail_mesh),
    vtail = make_vortex_rings(vtail_mesh)
);

## Case
fs  = Freestream(
    alpha = 4.0, 
    beta  = 0.0, 
    omega = [0., 0., 0.]
);

ref = References(
    speed     = 150., 
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing)
)

## Solve
@time sys = solve_case(
    aircraft, fs, ref;
    compressible     = true, # Compressibility correction flag
    print            = true, # Prints the results for only the aircraft
    print_components = true, # Prints the results for all components
);

## Compute forces, moments and velocities over each surface
ax_sys = Wind() # Axis systems: Geometry(), Stability(), Body()
@time CFs, CMs = surface_coefficients(sys; axes = ax_sys) # Coefficients
# Fs, Ms   = surface_dynamics(sys; axes = ax) # Forces and moments
# Fs       = surface_forces(sys; axes = ax) # Forces only
# vels     = surface_velocities(sys) # Velocities

## Aerodynamic coefficients
nf  = nearfield(sys)
ff  = farfield(sys)

nfs = nearfield_coefficients(sys)
ffs = farfield_coefficients(sys)
    
## Force/moment coefficients and derivatives
@time dvs = freestream_derivatives(sys; 
    axes = ax_sys,
    print = true,
    print_components = true,
    # farfield = true
);

## Viscous drag prediction

# Equivalent flat-plate skin friction estimation
CDv_wing  = profile_drag_coefficient(wing,  [0.8, 0.8], sys.reference)
CDv_htail = profile_drag_coefficient(htail, [0.6, 0.6], sys.reference)
CDv_vtail = profile_drag_coefficient(vtail, [0.6, 0.6], sys.reference)

CDv_plate = CDv_wing + CDv_htail + CDv_vtail

## Local dissipation form factor friction estimation
import LinearAlgebra: norm

edge_speeds = norm.(surface_velocities(sys)); # Inviscid speeds on the surfaces

# Drag coefficients
CDvd_wing  = profile_drag_coefficient(wing_mesh,  [0.8, 0.8], edge_speeds.wing,  sys.reference)
CDvd_htail = profile_drag_coefficient(htail_mesh, [0.6, 0.6], edge_speeds.htail, sys.reference)
CDvd_vtail = profile_drag_coefficient(vtail_mesh, [0.6, 0.6], edge_speeds.vtail, sys.reference)

CDv_diss = CDvd_wing + CDvd_htail + CDvd_vtail

## Viscous drag coefficient
CDv = CDv_diss

## Total force coefficients with empirical viscous drag prediction
CDi_nf, CY_nf, CL_nf, Cl, Cm, Cn = nf = nearfield(sys) 
CDi_ff, CY_ff, CL_ff = ff = farfield(sys)

nf_v = [ CDi_nf + CDv; CDv; nf ]
ff_v = [ CDi_ff + CDv; CDv; ff ]

print_coefficients(nf_v, ff_v)

## Plotting
#=========================================================#

using Plots
gr()

using LaTeXStrings
const LS = LaTeXString

## Coordinates
Plots.plot(
    aspect_ratio = 1,
    camera = (30, 30),
    zlim = span(wing) .* (-0.5, 0.5),
    size = (800, 600)
)
Plots.plot!(wing_mesh, label = "Wing")
Plots.plot!(htail_mesh, label = "Wing")
Plots.plot!(vtail_mesh, label = "Wing")
Plots.plot!(sys, wing_mesh, dist = 10, num_stream = 50, span = 10, color = :green)
Plots.plot!(sys, htail_mesh, dist = 2, num_stream = 50, span = 10, color = :green)
Plots.plot!(sys, vtail_mesh, dist = 2, num_stream = 50, span = 10, color = :green)

## Compute spanwise loads
wing_ll = spanwise_loading(wing_mesh, ref, CFs.wing,  sys.circulations.wing)
htail_ll = spanwise_loading(htail_mesh, ref, CFs.htail, sys.circulations.htail)
vtail_ll = spanwise_loading(vtail_mesh, ref, CFs.vtail, sys.circulations.vtail);

## Plot spanwise loadings
plot_CD = begin
    Plots.plot(label = :none, ylabel = "CDi")
    Plots.plot!(wing_ll[:,1], wing_ll[:,2], label = "Wing")
    Plots.plot!(htail_ll[:,1], htail_ll[:,2], label = "Horizontal Tail")
    Plots.plot!(vtail_ll[:,1], vtail_ll[:,2], label = "Vertical Tail")
end

plot_CY = begin
    Plots.plot(label = :none, ylabel = "CY")
    Plots.plot!(wing_ll[:,1], wing_ll[:,3], label = "Wing")
    Plots.plot!(htail_ll[:,1], htail_ll[:,3], label = "Horizontal Tail")
    Plots.plot!(vtail_ll[:,1], vtail_ll[:,3], label = "Vertical Tail")
end

plot_CL = begin
    Plots.plot(label = :none, ylabel = "CL")
    Plots.plot!(wing_ll[:,1], wing_ll[:,4], label = "Wing")
    Plots.plot!(htail_ll[:,1], htail_ll[:,4], label = "Horizontal Tail")
    Plots.plot!(vtail_ll[:,1], vtail_ll[:,4], label = "Vertical Tail")
    Plots.plot!(wing_ll[:,1], wing_ll[:,5], label = "Wing Normalized", xlabel = "y")
end

Plots.plot(plot_CD, plot_CY, plot_CL, size = (800, 700), layout = (3,1))

## VARIABLE ANALYSES
#=========================================================#

using Accessors
using Base.Iterators: product

## Speed variation
Vs = 1.0:10:300
res_Vs = combinedimsview(
    map(Vs) do V
        ref1 = @set ref.speed = V
        sys = solve_case(aircraft, fs, ref1, compressible = true)
        [ mach_number(ref1); farfield(sys); nearfield(sys) ]
    end, (1)
)

## Plot
Plots.plot(
    res_Vs[:,1], res_Vs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "M",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)

## Alpha variation
αs = -5:0.5:5
res_αs = combinedimsview(
    map(αs) do α
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref, compressible = true)
        [ α; farfield(sys); nearfield(sys) ]
    end, (1)
)

## Plot
Plots.plot(
    res_αs[:,1], res_αs[:,2:end],
    layout = (3,3), size = (900, 800),
    xlabel = "α",
    labels = ["CD_ff" "CY_ff" "CL_ff" "CD" "CY" "CL" "Cl" "Cm" "Cn"]
)


## (Speed, alpha) variation
res = combinedimsview(
    map(product(Vs, αs)) do (V, α)
        ref1 = @set ref.speed = V
        fst = @set fs.alpha = deg2rad(α)
        sys = solve_case(aircraft, fst, ref1, compressible = true)
        [ mach_number(ref1); α; farfield(sys); nearfield(sys) ]
    end
)

res_p = permutedims(res, (3,1,2))

## CDi
plt_CDi_ff = Plots.plot(camera = (75, 30), ylabel = L"\alpha", xlabel = "V", zlabel = L"C_{D_i}")
[ Plots.plot!(res_p[:,1,n], res_p[:,2,n], res_p[:,3,n], label = "", c = :black) for n in axes(res_p,3) ]

# CL
plt_CL_ff = Plots.plot(camera = (75,15), 
ylabel = L"\alpha", xlabel = "V", zlabel = L"C_{L}", 
label = "", c = :black)
[ Plots.plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,8,n], label = "", c = :black) for n in axes(res_p,3) ]

# Cm
plt_Cm_ff = Plots.plot(camera = (100,30), ylabel = L"\alpha", xlabel = "V", zlabel = L"C_m", label = "")
[ Plots.plot!(
    res_p[:,1,n], res_p[:,2,n], res_p[:,10,n], label = "", c = :black) for n in axes(res_p,3) ]

# Plot
p = Plots.plot(plt_CDi_ff, plt_CL_ff, plt_Cm_ff, layout = (1,3), size = (1300, 400))

## Save
savefig(p, "coeffs.png")

## Fancy plotting using Makie, if needed
#==================================================================#

include("vlm_aircraft_makie_plot.jl")

## Speed variation
f2 = Figure()
ax1 = Axis(f2[1,1], xlabel = L"M", ylabel = L"C_{D_{ff}}")
lines!(res_Vs[:,1], res_Vs[:,2]) # CD_ff
ax2 = Axis(f2[1,2], xlabel = L"M", ylabel = L"C_{Y_{ff}}")
lines!(res_Vs[:,1], res_Vs[:,3])
ax3 = Axis(f2[1,3], xlabel = L"M", ylabel = L"CL_{ff}")
lines!(res_Vs[:,1], res_Vs[:,4])
ax4 = Axis(f2[2,1], xlabel = L"M", ylabel = L"C_D")
lines!(res_Vs[:,1], res_Vs[:,5])
ax5 = Axis(f2[2,2], xlabel = L"M", ylabel = L"C_Y")
lines!(res_Vs[:,1], res_Vs[:,6])
ax6 = Axis(f2[2,3], xlabel = L"M", ylabel = L"C_L")
lines!(res_Vs[:,1], res_Vs[:,7])
ax7 = Axis(f2[3,1], xlabel = L"M", ylabel = L"C_\ell")
lines!(res_Vs[:,1], res_Vs[:,8])
ax8 = Axis(f2[3,2], xlabel = L"M", ylabel = L"C_m")
lines!(res_Vs[:,1], res_Vs[:,9])
ax9 = Axis(f2[3,3], xlabel = L"M", ylabel = L"C_n")
lines!(res_Vs[:,1], res_Vs[:,10])

f2

## Alpha variation
f3 = Figure()
ax1 = Axis(f3[1,1], xlabel = L"α", ylabel = L"C_{D_{ff}}")
lines!(res_αs[:,1], res_αs[:,2]) # CD_ff
ax2 = Axis(f3[1,2], xlabel = L"M", ylabel = L"C_{Y_{ff}}")
lines!(res_αs[:,1], res_αs[:,3])
ax3 = Axis(f3[1,3], xlabel = L"M", ylabel = L"CL_{ff}")
lines!(res_αs[:,1], res_αs[:,4])
ax4 = Axis(f3[2,1], xlabel = L"M", ylabel = L"C_D")
lines!(res_αs[:,1], res_αs[:,5])
ax5 = Axis(f3[2,2], xlabel = L"M", ylabel = L"C_Y")
lines!(res_αs[:,1], res_αs[:,6])
ax6 = Axis(f3[2,3], xlabel = L"M", ylabel = L"C_L")
lines!(res_αs[:,1], res_αs[:,7])
ax7 = Axis(f3[3,1], xlabel = L"M", ylabel = L"C_\ell")
lines!(res_αs[:,1], res_αs[:,8])
ax8 = Axis(f3[3,2], xlabel = L"M", ylabel = L"C_m")
lines!(res_αs[:,1], res_αs[:,9])
ax9 = Axis(f3[3,3], xlabel = L"M", ylabel = L"C_n")
lines!(res_αs[:,1], res_αs[:,10])

f3