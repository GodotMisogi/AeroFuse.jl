## Visualisation for `doublet_source_aircraft.jl` — aircraft surfaces coloured by pressure
## coefficient, the slender-body fuselage line, and streamlines through the coupled field.
## Needs a Makie backend (e.g. `using Pkg; Pkg.add("CairoMakie")`).
include("doublet_source_aircraft.jl")

using StaticArrays
using LinearAlgebra
import AeroFuse.PanelGeometry: p1, p2, p3, p4
using CairoMakie
CairoMakie.activate!()

cps = surface_coefficients(sys)

# RK4 streamline integration through the solved field
function streamline(seed, ds, n)
    pts = Vector{SVector{3,Float64}}(undef, n + 1); pts[1] = seed
    for k in 1:n
        r  = pts[k]
        k1 = normalize(field_velocity(sys, r));            k2 = normalize(field_velocity(sys, r + 0.5ds*k1))
        k3 = normalize(field_velocity(sys, r + 0.5ds*k2)); k4 = normalize(field_velocity(sys, r + ds*k3))
        pts[k+1] = r + ds/6 * (k1 + 2k2 + 2k3 + k4)
    end
    pts
end

seeds   = [ SVector(-1.6, y, z) for y in LinRange(-3.6, 3.6, 9) for z in (-0.25, 0.0, 0.25, 0.6) ]
streams = [ streamline(s, 0.04, 220) for s in seeds ]

fig   = Figure(size = (1500, 900))
scene = LScene(fig[1, 1]; show_axis = false)

local surf
for s in keys(sys.surfaces)
    p = sys.surfaces[s]; cp = cps[s]; nc, ns = size(p)
    verts = Point3f[]; faces = Vector{NTuple{3,Int}}(); vcol = Float32[]
    for j in 1:ns, i in 1:nc
        b = length(verts)
        push!(verts, Point3f(p1(p[i,j])), Point3f(p2(p[i,j])), Point3f(p3(p[i,j])), Point3f(p4(p[i,j])))
        push!(faces, (b+1, b+2, b+3), (b+1, b+3, b+4))
        append!(vcol, fill(clamp(Float32(cp[i,j]), -1.2f0, 1f0), 4))
    end
    fmat = reduce(vcat, [ [f[1] f[2] f[3]] for f in faces ])
    global surf = mesh!(scene, verts, fmat; color = vcol, colormap = :RdBu, colorrange = (-1.2, 1), shading = NoShading)
end
Colorbar(fig[1, 2], surf, label = L"Pressure coefficient $C_p$", height = Relative(0.6))

# Fuselage surface: the actual HyperEllipse body, built from its surface `coordinates` grid
# (indices [circumferential, section, xyz]) as a translucent grey mesh — see the `Plots`
# recipe `fuselage_plot(::HyperEllipseFuselage)` in src/Tools/plot_tools.jl. `fuse` is the
# original geometry from the included case; the solved `sys.fuselage` is only its axis line.
let n_secs = 20, n_circ = 24
    coo = coordinates(fuse, LinRange(0, 1, n_secs), n_circ)   # (n_circ, n_secs*3, 3)
    nci, nse = size(coo, 1), size(coo, 2)
    verts = Point3f[]; faces = Vector{NTuple{3,Int}}()
    for k in 1:nse, i in 1:nci
        push!(verts, Point3f(coo[i, k, 1], coo[i, k, 2], coo[i, k, 3]))
    end
    idx(i, k) = (k - 1) * nci + i                            # column-major into `verts`
    for k in 1:nse-1, i in 1:nci-1
        push!(faces, (idx(i, k), idx(i+1, k), idx(i+1, k+1)), (idx(i, k), idx(i+1, k+1), idx(i, k+1)))
    end
    fmat = reduce(vcat, [ [f[1] f[2] f[3]] for f in faces ])
    mesh!(scene, verts, fmat; color = (:grey70, 0.55), shading = NoShading)
end

# Fuselage slender-body singularity axis line
isnothing(sys.fuselage) || lines!(scene, [ Point3f(el.rc) for el in sys.fuselage ]; color = :black, linewidth = 3)

for pts in streams
    lines!(scene, Point3f.(pts); color = (:steelblue, 0.7), linewidth = 1.0)
end

update_cam!(scene.scene, Vec3f(6.0, -9.0, 5.0), Vec3f(2.0, 0.0, 0.0))
fig[0, :] = Label(fig, "Doublet-source aircraft + slender-body fuselage (α = 3°)", fontsize = 22)

save(joinpath(@__DIR__, "doublet_source_aircraft.png"), fig; px_per_unit = 2)
fig
