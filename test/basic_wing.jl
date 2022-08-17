##
using AeroFuse
using LinearAlgebra
using StaticArrays

## Geometry
wing = WingSection(
    area      = 19.0, # Projected area
    aspect    = 8.3, # Aspect ratio
    dihedral  = 0.0, # Dihedral angle (deg)
    sweep     = 10.0, # Sweep angle (deg)
    w_sweep   = 0.25, # Sweep location normalized to chords ∈ [0,1]
    taper     = 1.0, # Taper ratio
    root_foil = naca4(0,0,1,2), # Root airfoil
    tip_foil  = naca4(0,0,1,2), # Tip airfoil
    position  = [0.0, 0, 0.0], # Location (m)
    angle     = 0., # Angle of incidence (deg)
    axis      = [0., 1., 0.0], # Spanwise direction
    symmetry  = true # Symmetric about x-z plane
)

# Meshing
wing_mesh = WingMesh(wing, [10], 20)

surf_pts  = surface_coordinates(wing_mesh)
surf_pans = 
    # ComponentArray(
        # wing = 
        make_panels(surf_pts)
    # )
# surf_pans_view = @view permutedims(surf_pans)[:]

# Freestream velocity
α = 5.0
β = 0.0
Umag = 1.
fs = Freestream( α, β, zeros(3))
V∞ = Umag * velocity(fs)

wake_pans = [ wake_panel(surf_pans[:,i], 10., velocity(fs)) for i in axes(surf_pans, 2) ]

σs = -dot.(Ref(V∞), AeroFuse.PanelGeometry.collocation_point.(surf_pans))

##
@time system = solve_system(surf_pans, fs, 1.0e5);

##
Nc, Ns = size(surf_pans)
Nf = Nc * Ns
φs = permutedims(reshape(system.singularities[1:Nf], Ns, Nc))
# φs[:,1:end÷2] = φs[end:-1:1,1:end÷2]

φs

##
function surface_velocities(panels, φs, σs, fs)
    Q = velocity(fs)
    map(CartesianIndices(panels[1:end-1,1:end-1])) do ind
        i, j = ind.I

        # Streamwise-lateral-normal
        dφ_ds = -(φs[i+1,j] - φs[i,j]) / (distance(panels[i+1,j], panels[i,j]))
        dφ_dl = -(φs[i,j+1] - φs[i,j]) / (distance(panels[i,j+1], panels[i,j]))
        dφ_dn = σs[i,j]

        Q + SVector(dφ_ds, dφ_dl, dφ_dn)
    end
end

vels = surface_velocities(surf_pans, φs, σs, fs)

##
vs = AeroFuse.DoubletSource.surface_velocities(system)
@time cls, cps = surface_coefficients(system, projected_area(wing_mesh));
println("Σᵢ Clᵢ: $(sum(cls))")

## Plotting
using Plots

plotlyjs()

plt = Plots.plot(aspect_ratio = 1, zlim = (-0.5, 0.5) .* span(wing_mesh))

map(surf_pans) do pan
    xyz = combinedimsview(panel_coordinates(pan), (1))
    Plots.plot!(xyz[:,1], xyz[:,2], xyz[:,3], color = :grey, label = "")
end
map(wake_pans) do pan
    xyz = combinedimsview(panel_coordinates(pan), (1))
    Plots.plot!(xyz[:,1], xyz[:,2], xyz[:,3], color = :blue, label = "")
end
# Plots.scatter!(Tuple.(surf_pts)[:], markersize = 0.1)

plot!()

# savefig(plt)

##
