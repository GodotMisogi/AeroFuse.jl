
# ## Flight Dynamics

# ### Longitudinal Motion

# A standard three degrees-of-freedom rigid body model is used for performing flight dynamics analyses in the 2-dimensional longitudinal plane formulated as an initial-value problem. The coupled differential equations are shown in a canonical state-space representation subject to initial conditions $\mathbf x_0$, with the time-evolution of the state vector $\mathbf x$ driven by the forcing function $\mathbf f$.

# ```math
# \begin{aligned}
#     \frac{d\mathbf x}{dt} & = \mathbf f(\mathbf x, t) \\
#     \frac{d}{dt}
#     \begin{bmatrix}
#         u_b \\
#         w_b \\
#         Q \\
#         x_e \\
#         y_e \\
#         \Theta \\
#         \delta_e \\
#         m
#     \end{bmatrix} & =
#     \begin{bmatrix}
#         -Qw_b + (T - D(V_\infty)\cos\alpha + L(V_\infty, \delta_e)\sin\alpha - W\sin\Theta) / m \\
#         Qu_b + (  - D(V_\infty)\sin\alpha - L(V_\infty, \delta_e)\cos\alpha + W\cos\Theta) / m \\
#         (M_A(\alpha, \delta_e, \hat Q) - T\Delta_{zT}) / I_{yy} \\
#         u_b \cos\Theta - w_b \sin\Theta \\
#         u_b \sin\Theta + w_b \cos\Theta \\
#         Q \\ 
#         g(\mathbf x, t) \\
#         -c_T T
#     \end{bmatrix}
# \end{aligned}
# ```
# The fuel burn over time (viz. reduction of mass $m$) is computed using a specific fuel consumption value with a linear dependence on the thrust. A manual controller law for the elevator deflection angle $\delta_e$ can also be implemented by providing the function $g(\mathbf x, t)$.

# ### Full Space

# A standard six degrees-of-freedom rigid body model is used for performing flight dynamics analyses in 3 dimensions.