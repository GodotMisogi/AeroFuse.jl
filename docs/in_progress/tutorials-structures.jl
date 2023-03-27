
# # ## Euler-Bernoulli Beam Structural Analysis
# # The tubular beam's relevant properties, viz. the Young's (elastic) modulus $E$, shear modulus $G$, and torsional moment of inertia $J$ must be specified to define the stiffness matrix for its discretization with $n$ elements.

# ## Material properties
# E = 1. # Elastic modulus
# G = 1. # Shear modulus
# J = 2. # Torsional moment of inertia
# n = 2  # Number of sections

# ## Stiffness matrix
# K = bending_stiffness_matrix(
#     fill(E, 2), 
#     fill(G, 2),
#     fill(J, 2), 
#     :z         # Direction of deflection
# )

# # Fixed, hinged beam subjected to force and moment at the center.

# ## Stiffness matrix
# A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3

# ## Load vector
# b = [-1000, 1000, 0]    # F2, M2, M3

# ## Solution
# x = A \ b

# ## Forces
# F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]

# # Propped cantilever beam subjected to force at free end.

# ## Stiffness matrix
# A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2

# ## Load vector
# b = [10, 0, 0]

# ## Solution
# x = A \ b

# ## Forces
# F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]