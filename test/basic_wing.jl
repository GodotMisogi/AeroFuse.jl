##
using AeroMDAO
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(
    area      = 19.0, # Projected area
    aspect    = 8.3, # Aspect ratio
    dihedral  = 0.0, # Dihedral angle (deg)
    sweep     = 10.0, # Sweep angle (deg)
    w_sweep   = 0.25, # Sweep location normalized to chords ∈ [0,1]
    taper     = 0.4, # Taper ratio
    root_foil = naca4(0,0,1,2), # Root airfoil
    tip_foil  = naca4(0,0,1,2), # Tip airfoil
    position  = [0.0, 0, 0.0], # Location (m)
    angle     = 0., # Angle of incidence (deg)
    axis      = [0., 1., 0.0], # Spanwise direction
    symmetry  = true # Symmetric about x-z plane
)

# Meshing
wing_mesh = WingMesh(wing, [10], 10)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = ComponentArray(wing = make_panels(surf_pts))

# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 0.0
β = 0.0
Umag = 1.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

σs = -dot.(Ref(V∞), AeroMDAO.PanelGeometry.collocation_point.(surf_pans.wing))

##
system = solve_system(surf_pans, fs, 1.0e5);

##
Nc, Ns = size(surf_pans.wing)
Nf = Nc * Ns
φs = permutedims(reshape(system.singularities[1:Nf], Ns,  Nc))

##
function surface_velocities(panels :: AbstractMatrix{<:AbstractPanel3D}, φs :: AbstractMatrix{<:Real}, σs, fs)
    Q = velocity(fs)
    map(CartesianIndices(panels[1:end-1,1:end-1])) do ind
        i, j = ind.I

        # Streamwise-lateral-normal
        dφ_ds = -(φs[i+1,j] - φs[i,j]) / (distance(panels[i+1,j], panels[i,j]))
        dφ_dl = -(φs[i,j+1] - φs[i,j]) / (distance(panels[i,j+1], panels[i,j]))
        dφ_dn = σs[i,j]

        Q + SVector(dφ_ds, dφ_dl, dφ_dn)
    end
    # Δrx = @views distance.(panels[1:end-1,:], panels[2:end,:])
    # Δφx = @views φs[1:end-1,:] .- φs[2:end,:]

    # vxs = Δφx ./ Δrx

    # Δry = @views @. distance(panels[:,1:end-1]) - collocation_point(panels[:,2:end])
    # Δφx = @views φs[:,1:end-1] .- φs[:,2:end]



end

vels = surface_velocities(surf_pans.wing, φs, σs)
vxs, vys = surface_velocities(system)
@time cls, cps = surface_coefficients(system);
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots
plt_surfs = plot_panels(surf_pans)

plt = Plots.plot(aspect_ratio = 1, zlim = (-2.0, 3.0))

[ Plots.plot!(pan, color = :grey) for pan in plt_surfs ] 
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

##
