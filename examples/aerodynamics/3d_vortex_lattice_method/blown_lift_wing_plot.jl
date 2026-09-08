## Blown-lift visualization: spanwise loading + the deflected slipstream field
#
# Builds the same wing as `blown_lift_wing.jl` and renders two views with CairoMakie:
#   (a) spanwise lift distribution for the unpowered, axial-blown and powered-flap cases, and
#   (b) the slipstream velocity field in the propeller's chordwise plane, showing the jet tube
#       and its downward turning behind the flap.
using AeroFuse
using StaticArrays
using LinearAlgebra
using Accessors
using CairoMakie
CairoMakie.activate!()

## Geometry and case (see blown_lift_wing.jl)
#==========================================================================================#
wing = Wing(foils = fill(naca4(2,4,1,2),2), chords = [1.0,0.6], twists = [0.0,0.0],
            spans = [4.0], dihedrals = [5.0], sweeps = [5.0], symmetry = true)
wing_mesh  = WingMesh(wing, [12], 6)
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

# Mirrored tractor props at y = ±1.2
props = [PropellerDisk(center = SVector(-1.0, 1.2,0.0), axis = SVector(1.0,0.0,0.0), radius = 1.5, thrust = 8000.0, torque = 120.0, sense =  1.0),
         PropellerDisk(center = SVector(-1.0,-1.2,0.0), axis = SVector(1.0,0.0,0.0), radius = 1.5, thrust = 8000.0, torque = 120.0, sense = -1.0)]
turn_props = [ auto_turn(p, flapped, ref) for p in props ] # Jet turning auto-derived from the flap

## Analyses
#==========================================================================================#
sys_clean   = VortexLatticeSystem(ComponentVector(wing = wing_rings), fs, ref)
sys_flap    = VortexLatticeSystem(ComponentVector(wing = flapped),    fs, ref)
sys_powered = VortexLatticeSystem(ComponentVector(wing = flapped),    fs, ref; slipstream = turn_props)

# Spanwise normalized lift loading (wind axes, column 5 = CL_norm from 2Γ/ρVc)
function spanload(sys)
    CFs, _ = surface_coefficients(sys; axes = Wind())
    L = spanwise_loading(wing_mesh, ref, CFs.wing, sys.circulations.wing)
    return L[:, 1], L[:, 5] # ys, CL_norm
end
y_c, cl_c = spanload(sys_clean)
y_f, cl_f = spanload(sys_flap)
y_p, cl_p = spanload(sys_powered)

## Slipstream velocity field in the propeller's chordwise (x–z) plane at y = +1.2
#==========================================================================================#
prop_r = turn_props[1]
y_slice = prop_r.center[2]
xs = range(-1.6, 3.0; length = 30)
zs = range(-1.3, 1.3; length = 18)

# Plot the slipstream *increment* (not the total flow): its magnitude is small next to the
# freestream, so the total barely tilts, but the increment alone shows the tube crisply and its
# downward turning behind the flap. A uniform tube has near-constant speed, so colour the
# arrows by the local jet turning angle (0 ahead of the flap → full turn behind it), which is
# what the deflected-slipstream model actually varies. Only in-tube points carry an arrow.
pts_x = Float64[]; pts_z = Float64[]; vec_x = Float64[]; vec_z = Float64[]; turn_deg = Float64[]
for x in xs, z in zs
    vs = slipstream_velocity(SVector(x, y_slice, z), prop_r, ref) # Slipstream increment (geometry axes)
    norm(vs) < 0.5 && continue
    push!(pts_x, x); push!(pts_z, z)
    push!(vec_x, vs[1]); push!(vec_z, vs[3]); push!(turn_deg, atand(-vs[3], vs[1]))
end

# Developed tube radius (same model as slipstream_velocity), for the tube outline
Vs = sqrt(ref.speed^2 + 2prop_r.thrust / (ref.density * π * prop_r.radius^2))
Rs = prop_r.radius * sqrt((ref.speed + Vs) / (2Vs))
xc = prop_r.center[1]; zc = prop_r.center[3]

## Figure
#==========================================================================================#
fig = Figure(size = (1180, 500))

# (a) Spanwise loading
ax1 = Axis(fig[1,1], xlabel = L"y \;[\mathrm{m}]", ylabel = L"c_\ell \, c / c_{ref}",
           title = "Spanwise lift loading")
vspan!(ax1, prop_r.center[2] - Rs, prop_r.center[2] + Rs, color = (:gray, 0.12))
vspan!(ax1, -prop_r.center[2] - Rs, -prop_r.center[2] + Rs, color = (:gray, 0.12))
lines!(ax1, y_c, cl_c, label = "Clean, unpowered",       linewidth = 2.5)
lines!(ax1, y_f, cl_f, label = "Flap 20°, unpowered",    linewidth = 2.5)
lines!(ax1, y_p, cl_p, label = "Flap 20°, powered (turned jet)", linewidth = 2.5)
axislegend(ax1, position = :cb, framevisible = false)

# (b) Slipstream field with the turned jet
ax2 = Axis(fig[1,2], xlabel = L"x \;[\mathrm{m}]", ylabel = L"z \;[\mathrm{m}]",
           title = "Slipstream at y = $(round(y_slice, digits=1)) m", aspect = DataAspect())
arrows2d!(ax2, pts_x, pts_z, vec_x, vec_z;
          lengthscale = 0.009, color = turn_deg, colormap = :viridis)
lines!(ax2, [xc, xc], [zc - prop_r.radius, zc + prop_r.radius], color = :red, linewidth = 3) # Disk
lines!(ax2, [xc, 3.0], [zc + Rs, zc + Rs], color = (:red, 0.5), linestyle = :dash)            # Tube edges
lines!(ax2, [xc, 3.0], [zc - Rs, zc - Rs], color = (:red, 0.5), linestyle = :dash)
lines!(ax2, [0.0, 0.67, 1.0], [0.0, 0.0, -tan(δ_flap) * 0.33], color = :black, linewidth = 4) # Wing + flap
Colorbar(fig[1,3], limits = (minimum(turn_deg), maximum(turn_deg)), colormap = :viridis, label = "jet turning [°]")

## Save
outfile = joinpath(@__DIR__, "blown_lift_wing.png")
save(outfile, fig)
println("Saved figure to ", outfile)
println("CL: clean=$(round(nearfield(sys_clean).CZ,digits=3))  flap=$(round(nearfield(sys_flap).CZ,digits=3))  powered-flap=$(round(nearfield(sys_powered).CZ,digits=3))")
