using AeroFuse
using LinearAlgebra
using Accessors
using StaticArrays

function deflect_control_surfaces(rings, δ, h_l, is_control_surface)
    map(rings, is_control_surface) do ring, is_cs
        if is_cs
            n_l = cross(h_l, ring.normal) 
            new_normal = normalize(ring.normal + δ * n_l)
            @set ring.normal = new_normal
        else
            ring
        end
    end
end

wing = Wing(
    foils = fill(naca4(2, 4, 1, 2), 2),
    chords = [1.0, 0.6],
    spans = [4.0],
    symmetry = true
)
wing_mesh = WingMesh(wing, [12], 6)
wing_rings = make_vortex_rings(wing_mesh)

nc, ns = size(wing_rings)
half_ns = ns ÷ 2
is_aileron_right = zeros(Bool, size(wing_rings))
is_aileron_left = zeros(Bool, size(wing_rings))

is_aileron_left[end-1:end, 1:half_ns÷2] .= true
is_aileron_right[end-1:end, (half_ns + half_ns÷2 + 1):end] .= true

δ_aileron = deg2rad(10.0)
wing_rings = deflect_control_surfaces(wing_rings, δ_aileron, SVector(0.0, 1.0, 0.0), is_aileron_right)
wing_rings = deflect_control_surfaces(wing_rings, -δ_aileron, SVector(0.0, 1.0, 0.0), is_aileron_left)

aircraft = ComponentVector(wing = wing_rings)
fs = Freestream(alpha=3.0, beta=0.0, omega=[0.0, 0.0, 0.0])
ref = References(speed=50.0, density=1.225, viscosity=1.5e-5, area=projected_area(wing), span=span(wing), chord=mean_aerodynamic_chord(wing), location=mean_aerodynamic_center(wing))

sys = solve_case(aircraft, fs, ref)
nf = nearfield(sys)
println("Nearfield: ", nf)

using Plots
gr()

println("Generating plot with streamlines...")
p = Plots.plot(
    aspect_ratio = 1,
    camera = (30, 30),
    zlim = span(wing) .* (-0.5, 0.5),
    size = (800, 600),
    title = "Aileron Deflection Streamlines",
)

Plots.plot!(p, wing_mesh, label = "Wing")
Plots.plot!(p, sys, wing_mesh, dist = 5, num_stream = 50, span = 10, color = :green)

savefig(p, "vlm_streamlines.png")
println("Saved plot to vlm_streamlines.png")
