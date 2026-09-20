using AeroFuse
using LinearAlgebra
using Accessors
using StaticArrays

# Helper to modify VortexRings to include a normal deflection using `transform_normal`
# h_l is the hinge axis unit vector, g_l is the blending gain 
function deflect_control_surfaces(rings, δ, h_l, is_control_surface)
    map(rings, is_control_surface) do ring, is_cs
        if is_cs
            # Drela's small angle normal transformation: n1 ≈ n0 + δ * (h_l × n0)
            n_l = cross(h_l, ring.normal) 
            # Equivalently, transform_normal(panel, h_l, 1.0) would return this.
            # Using exact rotation for arbitrarily large deflections is better:
            # but we use the linear approximation from Drela for demonstration.
            new_normal = normalize(ring.normal + δ * n_l)
            @set ring.normal = new_normal
        else
            ring
        end
    end
end

## Wing
wing = Wing(
    foils     = fill(naca4(2, 4, 1, 2), 2),
    chords    = [1.0, 0.6],
    twists    = [0.0, 0.0],
    spans     = [4.0],
    dihedrals = [5.0],
    sweeps    = [5.0],
    symmetry  = true,
)

wing_mesh = WingMesh(wing, [12], 6)
wing_rings = make_vortex_rings(wing_mesh)

# Define boolean array for control surfaces: Ailerons on the outer half of the span, trailing 2 chords
is_aileron_right = zeros(Bool, size(wing_rings))
is_aileron_left = zeros(Bool, size(wing_rings))

# Note: symmetry generates the negative y side first, then the positive y side
# size(wing_rings) is (num_chord, num_span*2) = (6, 24)
nc, ns = size(wing_rings)
half_ns = ns ÷ 2

# Left aileron
is_aileron_left[end-1:end, 1:half_ns÷2] .= true
# Right aileron 
is_aileron_right[end-1:end, (half_ns + half_ns÷2 + 1):end] .= true

# Deflect ailerons: trailing edge down on right wing (δ > 0), trailing edge up on left wing (δ < 0)
hinge_axis_right = SVector(0.0, 1.0, 0.0) # approx spanwise
hinge_axis_left  = SVector(0.0, -1.0, 0.0) # approx spanwise, flip for symmetry

δ_aileron = deg2rad(10.0)

# Modify the normals using the linear transform
wing_rings = deflect_control_surfaces(wing_rings, δ_aileron, hinge_axis_right, is_aileron_right)
# Left goes opposite direction for roll
wing_rings = deflect_control_surfaces(wing_rings, -δ_aileron, hinge_axis_left, is_aileron_left)

aircraft = ComponentVector(wing = wing_rings)

## Case
fs = Freestream(alpha = 3.0, beta = 0.0, omega = [0.0, 0.0, 0.0])

ref = References(
    speed     = 50.0,
    density   = 1.225,
    viscosity = 1.5e-5,
    area      = projected_area(wing),
    span      = span(wing),
    chord     = mean_aerodynamic_chord(wing),
    location  = mean_aerodynamic_center(wing),
)

## Solve
sys = solve_case(
    aircraft, fs, ref;
    compressible     = true,
    print            = true,
    print_components = true,
)

ax = Wind()
CFs, CMs = surface_coefficients(sys; axes = ax)

# Expect nonzero rolling moment Cl due to antisymmetrical aileron deflection
nf = nearfield(sys)
ff = farfield(sys)
println("Nearfield (CD, CY, CL, Cl, Cm, Cn) = ", nf)
println("Rolling Moment (Cl): ", nf[4])
