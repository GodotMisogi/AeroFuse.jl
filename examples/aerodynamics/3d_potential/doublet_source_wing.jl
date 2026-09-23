## 3D doublet-source (Morino) panel method — NACA 0012 rectangular wing
using AeroFuse
using LinearAlgebra

# Geometry: symmetric NACA 0012 rectangular wing, aspect ratio 8
foil = naca4(0, 0, 1, 2)
wing = Wing(
    foils     = [foil, foil],
    chords    = [1.0, 1.0],
    spans     = [4.0],
    dihedrals = [0.0],
    sweeps    = [0.0],
    symmetry  = true,
)

@show span(wing), projected_area(wing), aspect_ratio(wing)

# Surface panel mesh (chordwise wraps around the airfoil, spanwise across the wing)
mesh      = WingMesh(wing, [8], 20)
surf_pans = surface_panels(mesh)
@show size(surf_pans)

# Freestream and solve
α  = 5.0
fs = Freestream(alpha = α, beta = 0.0)

sys = AeroFuse.solve_system(surf_pans, fs, 1e3)
println(sys)

# Aerodynamic coefficients
cls, cps = AeroFuse.surface_coefficients(sys, projected_area(wing))

println("CL   = ", round(sum(cls), digits = 4))
println("cp   ∈ [", round(minimum(cps), digits = 3), ", ", round(maximum(cps), digits = 3), "]")
println("lifting-line reference CL ≈ ",
        round(2π * deg2rad(α) / (1 + 2 / aspect_ratio(wing)), digits = 4))
