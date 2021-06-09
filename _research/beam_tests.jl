using AeroMDAO

## Deflection stiffness matrix
K = deflection_stiffness_matrix([1., 1.], [1., 1.], [2., 2.], :z)

## 1. Fixed hinged beam subjected to force and moment at the center
A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3 
b = [-1000, 1000, 0]    # F2, M2, M3

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
J = torsional_stiffness_matrix([1., 1., 1.], [1., 1., 1.], [2., 2., 2.])

## 1. ???
A = J[[1,2],[1,2]] # ψ1, ψ2 
b = [-1000, 1000]  # R2, R2

x = A \ b

M = J * [ x; zeros(2) ]

@test M ≈ [-1000., 1000., 0., 0.] atol = 1e-6