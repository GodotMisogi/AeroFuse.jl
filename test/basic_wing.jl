##
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(
    area      = 19.0, # Projected area
    aspect    = 8.3, # Aspect ratio
    dihedral  = 3.0, # Dihedral angle (deg)
    sweep     = 10.0, # Sweep angle (deg)
    w_sweep   = 0.25, # Sweep location normalized to chords ∈ [0,1]
    taper     = 0.4, # Taper ratio
    root_foil = naca4(0,0,1,2), # Root airfoil
    tip_foil  = naca4(0,0,1,2), # Tip airfoil
    position  = [4.0, 0, -0.5], # Location (m)
    angle     = 0., # Angle of incidence (deg)
    axis      = [0., 1., 0.0], # Spanwise direction
    symmetry  = true # Symmetric about x-z plane
)

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
system = solve_system(surf_pans, Umag, fs, 1.0e5);

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
	Δφx = φs[:,1:end-1] .- φs[:,2:end]
end


## Plotting
using Plots

plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.2, 1.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)