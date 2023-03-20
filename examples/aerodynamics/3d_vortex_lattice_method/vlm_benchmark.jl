##
using BenchmarkTools

## BYU Flow
using VortexLattice

function vlm_byu()
    # Simple Wing with Uniform Spacing
    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 6
    spacing_s = VortexLattice.Uniform()
    spacing_c = VortexLattice.Uniform()

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = VortexLattice.Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = VortexLattice.Freestream(Vinf, alpha, beta, Omega)

    # vortex rings with mirrored geometry
    mirror = true
    symmetric = false

    _, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    surfaces = [surface]

    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=VortexLattice.Stability())
    CDiff = far_field_drag(system)

    return CF, CM, CDiff, system
end 

##
@time CF, CM, CDiff, byu_sys = vlm_byu();

print_coefficients([CF; CM], [CDiff,"—","—"], "BYU")

## AeroFuse.jl
using AeroFuse

## Wing section setup
function vlm_aerofuse()
    # Wing with same settings
    wing = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [2.2, 1.8],
        twists    = [2.0, 2.0],
        spans     = [7.5],
        dihedrals = [0.],
        sweeps    = [3.0528],
        w_sweep   = 0., # Quarter-chord sweep
        symmetry  = true,
    )

    wing_mesh = WingMesh(wing, [24], 6, span_spacing = AeroFuse.Uniform());

    # Freestream conditions
    fs  = Freestream(
        alpha = 1.0, # deg
        beta  = 0.0, # deg
        omega = [0.,0.,0.]
    )

    # Reference values
    ref = References(
        speed     = 1., # m/s
        density   = 1.225, # kg/m³
        viscosity = 1.5e-5, # ???
        area      = 30.0, # m²
        span      = 15.0, # m
        chord     = 2.0, # m
        location  = [0.50, 0.0, 0.0] # m
    )

    ## Horseshoes
    ac_hs = ComponentVector(wing = make_horseshoes(wing_mesh))
    system = VortexLatticeSystem(ac_hs, fs, ref)

    ## Vortex rings
    # ac_vs =  ComponentVector(wing = make_vortex_rings(wing_mesh))
    # system = VortexLatticeSystem(ac_vs, fs, ref)

    return nearfield(system), farfield(system), system
end

##
@time nfs, ffs, sys = vlm_aerofuse();

print_coefficients(nfs, ffs, "AeroFuse")

## Results
# ───────────────────────────────────
#  BYU  Nearfield          Farfield  
# ───────────────────────────────────
#  CX   0.00246661  CDi   0.00247615
#  CY      0.0      CYff      —
#  CZ    0.244371    CL       —
#  Cl      0.0                —
#  Cm   -0.0208467            —
#  Cn      0.0                —
# ───────────────────────────────────

# ────────────────────────────────────────
#  AeroFuse  Nearfield          Farfield  
# ────────────────────────────────────────
#     CX     0.0023826   CDi   0.00249797
#     CY        -0.0     CYff     -0.0
#     CZ      0.245316    CL    0.245319
#     Cl        0.0                —
#     Cm     -0.0211525            —
#     Cn        0.0                —
# ────────────────────────────────────────

## Benchmarking
@benchmark vlm_byu()

# BenchmarkTools.Trial: 1938 samples with 1 evaluation.
#  Range (min … max):  2.348 ms … 53.391 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     2.393 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   2.578 ms ±  1.586 ms  ┊ GC (mean ± σ):  3.06% ± 6.92%

#   █▃▁                                                         
#   ███▇█▆▅▅▅▅▃▃▃▄▁▁▃▃▄▄▁▃▁▄▃▁▃▃▃▁▃▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▃▁▃▁▁▁▁▃▁▁▁▃ █
#   2.35 ms      Histogram: log(frequency) by time     8.82 ms <

#  Memory estimate: 941.44 KiB, allocs estimate: 8188.

##
@benchmark vlm_aerofuse()

# BenchmarkTools.Trial: 4353 samples with 1 evaluation.
#  Range (min … max):  1.037 ms …  14.134 ms  ┊ GC (min … max): 0.00% … 91.51%
#  Time  (median):     1.085 ms               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   1.147 ms ± 666.701 μs  ┊ GC (mean ± σ):  3.80% ±  6.19%

#      █▆▁                                                       
#   ▃▄▇███▇▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▁▂▁▂▂▂▂▁▁▂▁▁▂▂▁▁▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▂▁▁▂ ▃
#   1.04 ms         Histogram: frequency by time        1.67 ms <

#  Memory estimate: 587.05 KiB, allocs estimate: 1230.