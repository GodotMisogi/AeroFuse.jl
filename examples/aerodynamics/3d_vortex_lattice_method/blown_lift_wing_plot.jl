## Blown-lift visualization: spanwise loading + the deflected slipstream field
#
# Builds a flapped wing and renders three views with CairoMakie:
#   (a) circulation loading (2Γ/ρV∞c): an axial slipstream *reduces* the bound circulation
#       under the prop (lower effective incidence), while a turned jet *raises* it,
#   (b) force loading (local dynamic pressure): both rise under the prop, because the higher
#       slipstream dynamic pressure lifts the force even where the circulation fell, and
#   (c) the slipstream velocity field, showing the jet tube and its downward turning past the
#       flap.
# Localized mirrored props are used so the under-prop features are distinct.
using AeroFuse
using StaticArrays
using LinearAlgebra
using Accessors
using CairoMakie
CairoMakie.activate!()

## Geometry and case
#==========================================================================================#
wing = Wing(foils = fill(naca4(2,4,1,2),2), chords = [1.0,0.6], twists = [0.0,0.0],
            spans = [4.0], dihedrals = [5.0], sweeps = [5.0], symmetry = true)
wing_mesh  = WingMesh(wing, [30], 16)
wing_rings = make_vortex_rings(wing_mesh)

fs  = Freestream(alpha = 3.0, beta = 0.0, omega = [0.0,0.0,0.0])
ref = References(speed = 50.0, density = 1.225, viscosity = 1.5e-5, area = projected_area(wing),
                 span = span(wing), chord = mean_aerodynamic_chord(wing), location = mean_aerodynamic_center(wing))

# Trailing-edge flaps over the full span, deflected 20° (Drela small-angle normal update)
δ_flap  = deg2rad(20.0)
nc, ns  = size(wing_rings)
is_flap = zeros(Bool, nc, ns); is_flap[end-1:end, :] .= true
flapped = map(wing_rings, is_flap) do ring, is_cs
    is_cs ? (@set ring.normal = normalize(ring.normal + δ_flap * cross(SVector(0.0,1.0,0.0), ring.normal))) : ring
end

# Localized mirrored tractor props at y = ±1.7
props        = [PropellerDisk(center = SVector(-1.0, 1.7,0.0), axis = SVector(1.0,0.0,0.0), radius = 1.0, thrust = 8000.0, torque = 120.0, sense =  1.0),
                PropellerDisk(center = SVector(-1.0,-1.7,0.0), axis = SVector(1.0,0.0,0.0), radius = 1.0, thrust = 8000.0, torque = 120.0, sense = -1.0)]
axial_props  = props                                         # No turning
turned_props = [ auto_turn(p, flapped, ref) for p in props ] # Jet turning auto-derived from the flap

## Analyses (all on the flapped wing)
#==========================================================================================#
sys_flap   = VortexLatticeSystem(ComponentVector(wing = flapped), fs, ref)
sys_axial  = VortexLatticeSystem(ComponentVector(wing = flapped), fs, ref; slipstream = axial_props)
sys_turned = VortexLatticeSystem(ComponentVector(wing = flapped), fs, ref; slipstream = turned_props)

# Spanwise loading: column 1 = y, column 4 = force-based CL (includes local q, wind axes),
# column 5 = circulation-based CL_norm = 2Γ/ρV∞c (no local-q term).
function spanload(sys)
    CFs, _ = surface_coefficients(sys; axes = Wind())
    L = spanwise_loading(wing_mesh, ref, CFs.wing, sys.circulations.wing)
    return L[:, 1], L[:, 5], L[:, 4] # ys, circulation loading, force loading
end
y, g_flap, f_flap     = spanload(sys_flap)
_, g_axial,  f_axial  = spanload(sys_axial)
_, g_turned, f_turned = spanload(sys_turned)

## Slipstream velocity field in the right prop's chordwise (x–z) plane
#==========================================================================================#
prop_r  = turned_props[1]
y_slice = prop_r.center[2]
xs = range(-1.6, 3.0; length = 30)
zs = range(-1.3, 1.3; length = 18)

# Plot the slipstream *increment* coloured by local jet turning angle (0 ahead of the flap →
# full turn behind it); only in-tube points carry an arrow.
pts_x = Float64[]; pts_z = Float64[]; vec_x = Float64[]; vec_z = Float64[]; turn_deg = Float64[]
for x in xs, z in zs
    vs = slipstream_velocity(SVector(x, y_slice, z), prop_r, ref)
    norm(vs) < 0.5 && continue
    push!(pts_x, x); push!(pts_z, z)
    push!(vec_x, vs[1]); push!(vec_z, vs[3]); push!(turn_deg, atand(-vs[3], vs[1]))
end

Vs = sqrt(ref.speed^2 + 2prop_r.thrust / (ref.density * π * prop_r.radius^2))
Rs = prop_r.radius * sqrt((ref.speed + Vs) / (2Vs))
xc = prop_r.center[1]; zc = prop_r.center[3]; yp = prop_r.center[2]

## Figure
#==========================================================================================#
fig = Figure(size = (1350, 640))

shade!(ax) = for s in (1, -1) # Mark both prop tubes
    vspan!(ax, s*yp - Rs, s*yp + Rs, color = (:gray, 0.12))
end
c_flap, c_axial, c_turned = :steelblue, :darkorange, :seagreen

# (a) Circulation loading
ax1 = Axis(fig[1,1], ylabel = L"2\Gamma / \rho V_\infty c", title = "Circulation loading")
shade!(ax1)
lines!(ax1, y, g_flap,   color = c_flap,   linewidth = 2.5, label = "Flap 20°, unpowered")
lines!(ax1, y, g_axial,  color = c_axial,  linewidth = 2.5, label = "+ axial jet")
lines!(ax1, y, g_turned, color = c_turned, linewidth = 2.5, label = "+ turned jet")
axislegend(ax1, position = :cb, framevisible = false, labelsize = 11)

# (b) Force loading
ax2 = Axis(fig[2,1], xlabel = L"y \;[\mathrm{m}]", ylabel = L"c_\ell\, c / c_{ref}\;\; (\mathrm{local}\ q)",
           title = "Force loading")
shade!(ax2)
lines!(ax2, y, f_flap,   color = c_flap,   linewidth = 2.5)
lines!(ax2, y, f_axial,  color = c_axial,  linewidth = 2.5)
lines!(ax2, y, f_turned, color = c_turned, linewidth = 2.5)
linkxaxes!(ax1, ax2); hidexdecorations!(ax1, grid = false)

# (c) Slipstream field with the turned jet
ax3 = Axis(fig[1:2,2], xlabel = L"x \;[\mathrm{m}]", ylabel = L"z \;[\mathrm{m}]",
           title = "Slipstream at y = $(round(y_slice, digits=1)) m", aspect = DataAspect())
arrows2d!(ax3, pts_x, pts_z, vec_x, vec_z; lengthscale = 0.005, color = turn_deg, colormap = :viridis)
lines!(ax3, [xc, xc], [zc - prop_r.radius, zc + prop_r.radius], color = :red, linewidth = 3)  # Disk
lines!(ax3, [xc, 3.0], [zc + Rs, zc + Rs], color = (:red, 0.5), linestyle = :dash)             # Tube edges
lines!(ax3, [xc, 3.0], [zc - Rs, zc - Rs], color = (:red, 0.5), linestyle = :dash)
lines!(ax3, [0.0, 0.67, 1.0], [0.0, 0.0, -tan(δ_flap) * 0.33], color = :black, linewidth = 4) # Wing + flap
Colorbar(fig[1:2,3], limits = (minimum(turn_deg), maximum(turn_deg)), colormap = :viridis, label = "jet turning [°]")

colsize!(fig.layout, 1, Relative(0.42))

## Save
outfile = joinpath(@__DIR__, "blown_lift_wing.png")
save(outfile, fig)
println("Saved figure to ", outfile)
println("CL: flap=$(round(nearfield(sys_flap).CZ,digits=3))  +axial=$(round(nearfield(sys_axial).CZ,digits=3))  +turned=$(round(nearfield(sys_turned).CZ,digits=3))")
