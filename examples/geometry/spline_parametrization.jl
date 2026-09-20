# %%
using AeroFuse
using Plots

## Wing section setup
wing = Wing(
    foils=fill(naca4((2, 4, 1, 2)), 3),  # Airfoils, type Foil
    chords=[2.0, 1.6, 0.8],            # Chord lengths (m)
    twists=[0.0, 0.0, 0.0],            # Twist angles (deg)
    spans=[3.0, 1.2],                  # Span lengths (m)
    dihedrals=[5.0, 5.0],                   # Dihedral angles (deg)
    sweeps=[20.0, 40.0],                  # Sweep angles (deg)
    # controls  = [Flap(10, 0.5), Aileron(5, 0.5)],
    sweep_ratio=0.25,                       # Chord length fraction of sweep location
    symmetry=true,                       # Symmetry in x-z plane
    # flip      = true                        # Reflection about x-z plane
)

##
x_w, y_w, z_w = wing_mac = mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Meshing and assembly
wing_mesh = WingMesh(
    wing,
    [24, 6],
    6,
    #  span_spacing = Cosine()
)

# %% Get the normalized spanwise coordinates of the mesh
function normalize_span_coordinates(xyzs)
    # Get the spanwise coordinates of the mesh
    ys = xyzs[1, :, 2]

    # Normalize the spanwise coordinates to [0, 1]
    ys_normalized = @. (ys[1:end] - ys[1]) / (ys[end] - ys[1])

    # xs = xyzs[:, :, 1] # Get the chordwise coordinates
    cs = xyzs[end, 1:end, 1] - xyzs[1, 1:end, 1]
    cs_normalized = cs ./ maximum(cs) # Normalize the chordwise coordinates

    return ys_normalized, cs_normalized
end

xyzs = combinedims(camber_coordinates(wing_mesh), (1, 2))
ys, cs = normalize_span_coordinates(xyzs)

# %% Define spline
using BSplineKit

n_cp = 10 # Number of control points
y_cp = LinRange(0, 1, n_cp) # Normalized spanwise coordinates
c_cp = ones(n_cp) # Control points for chord length scales
t_cp = zeros(n_cp) # Control points for twist angles

itp = BSplineKit.interpolate(y_cp, c_cp, BSplineOrder(3)) # Define a spline
# itp.spline.coefs[3:end-2] .= 0.2
itp.spline.coefs[end÷2 .+ (-2:2)] .= 5.0 # Modify the spline coefficients

chord_scale = itp.(ys)

# %%
plot(
    y_cp,
    c_cp,
    label="Original",
    ylabel="Chord Length Scale",
    xlabel="Normalized spanwise coordinate",
)
scatter!(
    ys,
    chord_scale,
    label="Modified (Interpolated)",
)

#%% Modify the mesh
using CoordinateTransformations, Rotations, LinearAlgebra, StaticArrays

function scale_mesh_chords!(new_xyzs, xyzs, chord_scales)
    xs = xyzs[:, :, 1] # This should be a copy
    @views le, te = xs[1, :], xs[end, :]
    ref_axis = (1 - 0.25) * le + 0.25 * te # Reference axis for scaling

    @. new_xyzs[:, :, 1] = (xs - ref_axis') * chord_scales' + ref_axis'

    return nothing
end

nc, ns, nd = size(xyzs) # Number of chordwise, spanwise, and coordinate dimension

scale_x(c_scale) = @SMatrix [ c_scale 0 0; 0 1 0; 0 0 1 ] # Scale function for x-coordinates

affs = AffineMap.(
    scale_x.(chord_scale), # Scale factor
    Ref(zeros(3)), # Translation vector
)

using Tullio
@tullio new_xyzs[i, j, k] = affs[j](xyzs[i, j, k]) # Apply the affine transformation
# new_xyzs = mapslices() ????

# %% Create a copy of the mesh coordinates
new_xyzs = copy(xyzs)
scale_mesh_chords!(new_xyzs, xyzs, chord_scale)

# %%
using ForwardDiff
using SparseArrays
using ComponentArrays

R = similar(xyzs)
x_in = ComponentArray(
    xyzs=xyzs,
    chord_scale=chord_scale
) # Input for the Jacobian
R_x = spzeros(length(R), length(x_in)); # For the Jacobian


# %%
ForwardDiff.jacobian!(R_x, (R, x) -> scale_mesh_chords!(R, x.xyzs, x.chord_scale), R, x_in);

# %%
using Plots
# plotlyjs()
gr()

plot(wing, 
    camera=(90, 90), 
    # zlim=(-0.5, 0.5) .* span(wing), 
    aspect_ratio=1
)

scatter!(
    new_xyzs[:, :, 1],
    new_xyzs[:, :, 2],
    new_xyzs[:, :, 3],
    label="",
    color=:cornflowerblue,
    alpha=0.5,
    ms=1,
)

# %% Aerodynamic analysis
aircraft = ComponentVector(wing = make_vortex_rings(splitdimsview(new_xyzs, (1, 2))))

# Freestream conditions
fs = Freestream(
    alpha = 2.0, # deg
    beta = 0.0, # deg
    omega = [0.0, 0.0, 0.0],
)

# Reference values
ref = References(
    speed = 150, # m/s
    density = 1.225, # kg/m³
    viscosity = 1.5e-5, # ???
    area = projected_area(wing), # m²
    span = span(wing), # m
    chord = mean_aerodynamic_chord(wing), # m
    location = mean_aerodynamic_center(wing), # m
)

## Solve system
@time system = solve_case(
        aircraft, fs, ref;
        compressible = true, # Compressibility correction option
        print            = true, # Prints the results for only the aircraft
        # print_components = true, # Prints the results for all components
    );