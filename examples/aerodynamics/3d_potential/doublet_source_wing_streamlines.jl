## Streamlines over a NACA 0012 wing from the 3D doublet-source (Morino) panel method.
## Companion to `doublet_source_wing.jl`. Needs a Makie backend, e.g.
##   julia> using Pkg; Pkg.add("CairoMakie")
using AeroFuse
using LinearAlgebra
using StaticArrays
import AeroFuse.PanelGeometry: p1, p2, p3, p4, collocation_point
import AeroFuse.Laplace: velocity
import AeroFuse.VortexLattice: bound_leg_velocity

using CairoMakie
CairoMakie.activate!()

# ------------------------------------------------------------------ solve
foil = naca4(0, 0, 1, 2)
wing = Wing(
    foils    = [foil, foil],
    chords   = [1.0, 1.0],
    spans    = [4.0],
    symmetry = true,
)

mesh      = WingMesh(wing, [12], 24)
surf_pans = surface_panels(mesh)
nc, ns    = size(surf_pans)
npanf     = nc * ns

fs  = Freestream(alpha = 5.0)
V∞  = velocity(fs)
sys = AeroFuse.solve_system(surf_pans, fs, 1e3)

cls, cps = AeroFuse.surface_coefficients(sys, projected_area(wing))
@show sum(cls)

# doublet strengths aligned with the (chord, span) panel grid, and the wake strengths
μ     = sys.singularities
μsurf = permutedims(reshape(μ[1:npanf], ns, nc))
μwake = μ[npanf+1:end]
wakes = sys.wake_panels

# ------------------------------------------------------- field velocity
# A constant-strength doublet panel is equivalent to a vortex ring of circulation
# Γ = μ around its four edges, so the induced velocity is a Biot–Savart sum.
function ring_velocity(pan, r)
    a1, a2, a3, a4 = p1(pan), p2(pan), p3(pan), p4(pan)
    bound_leg_velocity(r - a1, r - a2, 1.0) +
    bound_leg_velocity(r - a2, r - a3, 1.0) +
    bound_leg_velocity(r - a3, r - a4, 1.0) +
    bound_leg_velocity(r - a4, r - a1, 1.0)
end

function field_velocity(r)
    v = V∞
    @inbounds for j in 1:ns, i in 1:nc
        v += μsurf[i, j] * ring_velocity(surf_pans[i, j], r)
    end
    @inbounds for w in 1:ns
        v += μwake[w] * ring_velocity(wakes[w], r)
    end
    v
end

# RK4 streamline integration (arc-length parametrised)
function streamline(seed, ds, nsteps)
    pts = Vector{SVector{3,Float64}}(undef, nsteps + 1)
    pts[1] = seed
    for k in 1:nsteps
        r  = pts[k]
        k1 = normalize(field_velocity(r))
        k2 = normalize(field_velocity(r + 0.5ds*k1))
        k3 = normalize(field_velocity(r + 0.5ds*k2))
        k4 = normalize(field_velocity(r + ds*k3))
        pts[k+1] = r + ds/6 * (k1 + 2k2 + 2k3 + k4)
    end
    pts
end

# Seed a spanwise rake just upstream of the leading edge, clustered in z around the
# airfoil so the streamlines split over the upper and lower surfaces.
xs0   = -0.5
zs0   = LinRange(-0.3, 0.3, 17)
ys0   = LinRange(-3.5, 3.5, 6)
seeds = [ SVector(xs0, y, z) for y in ys0 for z in zs0 ]
streams = [ streamline(s, 0.02, 260) for s in seeds ]

# --------------------------------------------------------------- plot
fig   = Figure(size = (1400, 900))
scene = LScene(fig[1, 1]; show_axis = false)

# Wing surface coloured by pressure coefficient
verts  = Point3f[]
faces  = Vector{NTuple{3,Int}}()
vcols  = Float32[]
for j in 1:ns, i in 1:nc
    b = length(verts)
    push!(verts, Point3f(p1(surf_pans[i,j])), Point3f(p2(surf_pans[i,j])),
                 Point3f(p3(surf_pans[i,j])), Point3f(p4(surf_pans[i,j])))
    push!(faces, (b+1, b+2, b+3), (b+1, b+3, b+4))
    append!(vcols, fill(Float32(cps[i,j]), 4))
end
fmat = reduce(vcat, [ [f[1] f[2] f[3]] for f in faces ])

surf = mesh!(scene, verts, fmat; color = vcols, colormap = :RdBu, colorrange = (-2, 1),
             shading = NoShading)
Colorbar(fig[1, 2], surf, label = L"Pressure coefficient $C_p$", height = Relative(0.6))

# Streamlines (single subtle colour so the C_p surface stays the focus)
for pts in streams
    lines!(scene, Point3f.(pts); color = (:steelblue, 0.75), linewidth = 1.0)
end

# 3/4 front-above view (span running left-to-right, flow front-to-back)
update_cam!(scene.scene, Vec3f(4.0, -7.5, 4.5), Vec3f(0.5, 0.0, -0.1))

fig[0, :] = Label(fig, "NACA 0012 wing — doublet-source streamlines (α = 5°)", fontsize = 22)

save(joinpath(@__DIR__, "doublet_source_wing_streamlines.png"), fig; px_per_unit = 2)
fig
