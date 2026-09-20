## NOTE: This is called at the end of vlm_body_panel_aircraft.jl. It won't work by itself —
## it expects `sys`, `ref`, the surface meshes, and the half-wing geometries from that scope.
#
# Renders the coupled aircraft as a 3-D scene of surface-pressure contours, mirroring
# `vlm_aircraft_makie_plot.jl`:
#   * Fuselage skin (source panels)     → true one-sided surface pressure coefficient C_p, one
#     flat colour per panel (a constant-strength source panel has a piecewise-constant C_p).
#   * Lifting surfaces (vortex lattice)  → the camber lattice is a single zero-thickness sheet,
#     so it has no upper/lower skin; instead each panel carries a net sectional loading
#     ΔC_p = (F·n̂) S_ref / A. To read it like a two-sided pressure, the sheet is drawn twice,
#     offset ±δ along its normal, with the upper face coloured −ΔC_p/2 (suction) and the lower
#     face +ΔC_p/2 (pressure) — i.e. the same load split to opposite signs on opposing sides.
# Blue = suction (C_p < 0), red = pressure (C_p > 0) throughout, on a single shared diverging
# scale. Uses GLMakie: its real depth buffer resolves the two offset camber faces correctly from
# any orbit angle (top shows suction, bottom shows pressure), which CairoMakie's painter-order
# compositing cannot.
using GLMakie
GLMakie.activate!()

using LinearAlgebra: norm, dot, normalize
using StaticArrays: SVector

set_theme!()

## Lifting-surface sectional loading (signed ΔC_p per panel)
#=========================================================#
# Geometry-axis force coefficients, so the per-panel force dots cleanly with the mesh normals.
CFs, CMs = surface_coefficients(sys; axes = Geometry())
Sref = ref.area

# Net loading coefficient ΔC_p = (F·n̂) S_ref / A on every panel of a surface block. The stored
# `normal_vector` of a camber panel is the raw cross-product (‖·‖ = 2A), so it is normalised.
function delta_cp(CF_block, mesh)
    nrm = normal_vector.(camber_panels(mesh))
    A   = panel_area.(camber_panels(mesh))
    return map((f, n, a) -> dot(f, normalize(n)) * Sref / a, CF_block, nrm, A)
end

# Vertex normals consistent with the panel normals (extrapolated component-wise, then normalised),
# used to offset the two faces of the sheet.
function vertex_normals(mesh)
    nrm = normal_vector.(camber_panels(mesh))
    nx = extrapolate_point_mesh(getindex.(nrm, 1))
    ny = extrapolate_point_mesh(getindex.(nrm, 2))
    nz = extrapolate_point_mesh(getindex.(nrm, 3))
    return map((a, b, c) -> normalize(SVector(a, b, c)), nx, ny, nz)
end

# The `wing` block is hcat(port, starboard); split the loading back onto the two half meshes.
n_port = size(camber_coordinates(wing_l_mesh), 2) - 1
dcp = (
    wing_l = delta_cp(CFs.wing[:, 1:n_port],     wing_l_mesh),
    wing_r = delta_cp(CFs.wing[:, n_port+1:end], wing_r_mesh),
    htail  = delta_cp(CFs.htail, htail_mesh),
    vtail  = delta_cp(CFs.vtail, vtail_mesh),
)
lift_max = maximum(abs, vcat(vec.(values(dcp))...)) / 2   # max |±ΔC_p/2| on the lifting surfaces

## Body skin pressures (true one-sided C_p, flat per source panel)
#=========================================================#
# Each source panel → two triangles with its own duplicated corners, so its constant C_p renders
# as a flat facet (no interpolation across panel edges).
function body_surface_mesh(panels, cps)
    npan  = length(panels)
    pts   = Vector{SVector{3,Float64}}(undef, 4npan)
    cols  = Vector{Float64}(undef, 4npan)
    faces = Matrix{Int}(undef, 2npan, 3)
    for (k, (pan, cp)) in enumerate(zip(panels, cps))
        i0 = 4(k - 1)
        pts[i0+1], pts[i0+2], pts[i0+3], pts[i0+4] = pan.p1, pan.p2, pan.p3, pan.p4
        cols[i0+1:i0+4] .= cp
        faces[2k-1, :] = [i0 + 1, i0 + 2, i0 + 3]
        faces[2k,   :] = [i0 + 1, i0 + 3, i0 + 4]
    end
    return pts, faces, cols
end

body_cps = body_pressure_coefficients(sys)
body_pts, body_faces, body_cols = body_surface_mesh(sys.elements.body, body_cps)
body_max = maximum(abs, body_cols)

# One shared symmetric C_p scale for the whole aircraft (body + lifting surfaces). The body's
# stagnation/suction dominates the range, so the camber loading reads paler — which is honest: it
# is a genuinely smaller pressure excursion. Narrow this (e.g. 0.6) to boost lifting-surface contrast.
cp_lim = max(body_max, lift_max)
crange = (-cp_lim, cp_lim)

## Streamlines
#=========================================================#
# Seed from the exposed half-wing leading edges (the reference gross wing sits at the origin).
span_points = 16
dz = 1e-3
seed = [ chop_leading_edge(wing_l, span_points); chop_leading_edge(wing_r, span_points) ] .+ Ref(SVector(0., 0., dz))
streams = streamlines(sys, seed, 5.0, 80)

## Figure
#=========================================================#
δ = 0.006   # half-separation (m) between the two rendered faces of each camber sheet

# Draw a camber sheet as two faces offset ±δ along its normal: +n̂ side gets −ΔC_p/2 (suction),
# −n̂ side +ΔC_p/2 (pressure). GLMakie's depth buffer resolves which face is visible per pixel, so
# both sides read correctly from any orbit angle; the small offset only keeps the two coincident
# sheets from z-fighting.
function plot_two_sided!(scene, mesh, dcp_panels)
    V      = camber_coordinates(mesh)
    nv     = vertex_normals(mesh)
    cpv    = extrapolate_point_mesh(dcp_panels)
    connec = triangle_connectivities(LinearIndices(V))
    poly!(scene, vec(V .+ δ .* nv), connec, color = vec(-cpv ./ 2), colormap = :coolwarm, colorrange = crange)
    poly!(scene, vec(V .- δ .* nv), connec, color = vec( cpv ./ 2), colormap = :coolwarm, colorrange = crange)
end

fig = Figure(size = (1280, 720))
# Top-oblique view so the wing upper surfaces (suction, blue) are visible; a lower angle would
# instead show the lower surfaces (pressure, red). Axis decorations hidden for a clean scene.
ax  = Axis3(fig[1, 1], aspect = :data, elevation = 0.52, azimuth = -0.7, protrusions = 0)
hidedecorations!(ax); hidespines!(ax)

# Lifting surfaces (two-sided ΔC_p split)
plot_two_sided!(ax, wing_l_mesh, dcp.wing_l)
plot_two_sided!(ax, wing_r_mesh, dcp.wing_r)
plot_two_sided!(ax, htail_mesh,  dcp.htail)
plot_two_sided!(ax, vtail_mesh,  dcp.vtail)

# Fuselage skin (true C_p)
poly!(ax, body_pts, body_faces, color = body_cols, colormap = :coolwarm, colorrange = crange)

# Planform borders
lines!(ax, plot_planform(wing_l), color = :black)
lines!(ax, plot_planform(wing_r), color = :black)
lines!(ax, plot_planform(htail),  color = :black)
lines!(ax, plot_planform(vtail),  color = :black)

# Streamlines
for stream in eachcol(streams)
    lines!(ax, Point3f.(stream), color = (:grey40, 0.5), linewidth = 0.5)
end

# Single shared colourbar for the whole aircraft
Colorbar(fig[1, 2], colormap = :coolwarm, limits = crange, label = "Surface pressure coefficient, Cp")

fig[0, :] = Label(fig, "Coupled Body-Panel + Vortex-Lattice Surface Pressures", fontsize = 20)

fig

## Save figure
#=========================================================#
# save("plots/vlm_body_panel_pressures.png", fig, px_per_unit = 2)
