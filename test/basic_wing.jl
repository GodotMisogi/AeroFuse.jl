##
using Pkg
Pkg.activate(".")
using Revise
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(root_foil  = naca4(0,0,1,2),
                   tip_foil   = naca4(0,0,1,2),
                   root_chord = 1.0,
                   taper      = 1.0,
                   span       = 4.0,
                   dihedral   = 0.0,
                   sweep      = 0.0)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = make_panels(surf_pts)

# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 5.0
β = 0.0
Umag = 15.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

##
system = solve_system(surf_pans, Umag, fs, 1.0e5)

##
npancd, npansp = size(surf_pans)
npanf = npancd * npansp
φs = permutedims(reshape(system.singularities[1:npanf], npansp, npancd))

##
@views function panel_velocity(panels :: AbstractMatrix{<:AbstractPanel3D}, φs :: AbstractMatrix{<:Real})
    Δrx = distance.(panels[1:end-1,:], panels[2:end,:])
    Δφx = φs[1:end-1,:] .- φs[2:end,:]
    # vels = map((Δφx, Δrx, θ) -> Δφ / Δr + ifelse(sources, dot(u, θ), 0.), Δφs, Δrs, θs)
    vxs = Δφx ./ Δrx

    Δry = @. collocation_point(panels[:,1:end-1]) - collocation_point(panels[:,2:end])
    Δφy = φs[:,1:end-1] .- φs[:,2:end]

    vxs = Δφx ./ Δrx
end


## Plotting
using Plots

plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.2, 1.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)