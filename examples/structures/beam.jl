## Beam model
using AeroFuse

##
aluminum = Material(# Aluminum properties
    elastic_modulus = 85e9,  # Elastic modulus, N/m²
    shear_modulus = 25e9,  # Shear modulus, N/m²,
    yield_stress = 350e6, # Yield stress with factor of safety 2.5, N/m²,
    density = 1.6e3, # Density, kg/m³
)

##
n = 18 # Number of elements
Ls = 0:1/(n-1):1
beam = Beam(aluminum, Ls, (1e-1, 8e-2), (2e-2, 8e-3), true)
coords = [ LinRange(0., 0., n) LinRange(0., 1., n) LinRange(0., 0., n) ] # Coordinates


##
using Plots
gr()

thickness = 100.
normer(rs) = [ rs; rs[end] ] / maximum(rs)
r_norms = normer(beam.section.radius) * thickness

##
plot(
    coords[:,1], coords[:,2], coords[:,3], 
    color = :black,
    # m = (thickness, 0.8, :thermal, Plots.stroke(0)), 
    # linez = σ_ns,
    # c = :thermal,
    # cbar = true, 
    # label = ifelse(i == 1, "Deflected Beam Stresses", ""), 
    linestyle = :solid, 
    linewidth = beam.section.radius * 20
)
