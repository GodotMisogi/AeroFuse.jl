## Euler-Bernoulli beam tests
using AeroMDAO

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3 # Density, kg/m³
ν     = 0.3   # Poisson's ratio (UNUSED FOR NOW)
R     = 1e-2  # Outer radius, m
t     = 8e-3  # Thickness, m

# Create material and tubes
Ls       = [0.1, 0.1]
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

K = tube_stiffness_matrix(aluminum, tubes)

A = K[1:6,1:6]  # v2, φ2, φ3
b = [0., 0.1, 0.001, -0.1, 0.03, -0.3]    # F2, M2, M3

x = A \ b

## Forces
F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]

## 2. Propped cantilever beam with force at one end
A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2
b = [10, 0, 0]

x = A \ b

## Forces
F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]


## Torsional stiffness matrix
J = axial_stiffness_matrix([1., 1., 1.], [1., 1., 1.], [2., 2., 2.])

## 1. ???
A = J[[1,2],[1,2]] # ψ1, ψ2
b = [-1000, 1000]  # R2, R2

x = A \ b

M = J * [ x; zeros(2) ]

## Type tests
#==========================================================================================#

# Material - Steel: https://www.steelconstruction.info/Steel_material_properties
E     = 21e12  # Elastic modulus, N/m²
G     = 81e9   # Shear modulus, N/m²
σ_max = 275e6  # Yield stress, N/m²
ρ     = 8050.  # Density, kg/m³
ν     = 0.3    # Poisson's ratio
steel = Material(E, G, σ_max, ρ)

# Tube
L    = 1.   # Length, m
R    = 0.1  # Outer radius, m
t    = 0.01 # Thickness, m
tube = Tube(steel, L, R, t)

Rs = radii(tube)
A  = area(tube)
Iy = moment_of_inertia(tube)
Iz = moment_of_inertia(tube)
J  = polar_moment_of_inertia(tube)

## Stiffness matrix construction
n  = 10
x  = [ fill(E, n) fill(G, n) fill(A, n) fill(Iy, n) fill(Iz, n) fill(J, n) fill(L, n) ]

K1 = tube_stiffness_matrix(x)
K2 = tube_stiffness_matrix(steel, fill(tube, n))

## Plotting sparse matrices
using Plots
unicodeplots()

##
spy(K)

##
spy(K2)

## Differentiation tests
#==========================================================================================#

using ForwardDiff
using ReverseDiff

ForwardDiff.jacobian(tube_stiffness_matrix, x)