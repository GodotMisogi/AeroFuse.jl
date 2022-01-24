# # Theory

# The following theory described is meant to be minimal, with mainly the equations presented for reference when reading the code.

# ## Geometry

# ### Wing Parametrization
# 
# A **half-wing** is defined in terms of a nested trapezoid. A single trapezoid is called a **section**. A section consists of two foil profiles, and their associated chord lengths and twist angles. Between them is their span length with associated _leading-edge_ dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for foils, chords, and twists, and ``n - 1`` entries for spans, dihedrals, sweeps for some ``n \in \mathbb N``. The following illustration should help visualize the concept.
# 
# ![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

# ## Aerodynamics

# The aerodynamic analyses in AeroMDAO mainly utilize potential flow theory and solve problems using a _boundary element method_. This essentially is the following Laplace equation problem with the following Robin (?) boundary conditions:

# ```math
# \nabla^2 \phi = 0, \quad \mathbf V \equiv \nabla \phi \cdot \hat{\mathbf n} = 0, \quad \lim_{\mathbf r \to \infty} \phi(\mathbf r) \to 0
# ```

# Implementations for viscous-inviscid analyses are in progress.

# ### Doublet-Source Panel Method

# The doublet-source panel method predicts **inviscid, incompressible, irrotational, isentropic** external flow over surfaces in 2 dimensions. 

# Source and doublet singularities are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the corresponding free-field Green function for each singularity. 

# ```math

# ```


# As the system is linear, the velocities are added together to obtain the total induced velocity at a point $\mathbf r$.

# ```math

# ```

# ### Vortex Lattice Method

# The vortex lattice method predicts **inviscid, incompressible, irrotational, isentropic** external flow over "thin" surfaces in 3 dimensions. 

# Vortex filaments are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the Biot-Savart integral for a vortex line of length $\ell$ with a constant circulation strength $\Gamma$.

# ```math
# \mathbf V(\mathbf r, \mathbf r') = \frac{\Gamma}{4\pi} \int_0^\ell \frac{d\boldsymbol\ell' \times (\mathbf r - \mathbf r')}{|\mathbf r - \mathbf r'|^3}
# ```

# The vortices can be set up in various configurations consisting of bound or semi-infinite filaments, commonly in the form of _horseshoes_ or _vortex rings_.

# > 1. Horseshoe elements
# > 
# > 2. Vortex rings
# > 

# A quasi-steady freestream condition with velocity $\mathbf U$ and rotation $\boldsymbol\Omega$ (in the body's frame) is imposed for the external flow. The induced velocity at a point is given by:

# ```math
# \mathbf V_{\infty}(\mathbf r) = - (\mathbf U + \boldsymbol\Omega \times \mathbf r)
# ```

# As the system is linear, the velocities are added together to obtain the total induced velocity at a point $\mathbf r$.

# ```math
# \mathbf V(\mathbf r) = \sum_i \frac{\Gamma_i}{4\pi} \int_0^{\ell_i} \frac{d\boldsymbol\ell_i' \times (\mathbf r - \mathbf r_i')}{|\mathbf r - \mathbf r_i'|^3} + \mathbf V_\infty(\mathbf r)
# ```

# Imposing the Neumann boundary condition $\mathbf V \cdot \hat{\mathbf n} = 0$ defines the problem. The construction, in essence, fundamentally results in the following linear system to be solved:

# ```math
# \mathbf A \boldsymbol\Gamma = -V_{\infty} \cdot [\hat{\mathbf n}_i ]_{i = 1, \ldots, N}
# ```

# ## Structures

# The structural analyses in AeroMDAO utilize _linear finite-element methods_. 

# Particularly, a $2$-dimensional beam element model has been implemented following the standard formulation using cubic Hermite shape functions based on **Euler-Bernoulli beam theory**. These are embedded into a $3$-dimensional local coordinate system in the vortex lattice method without loss of generality.

# The linear system consisting of the stiffness matrix $\mathbf K$ and load vector $\mathbf f$ are solved to obtain the displacement vector $\boldsymbol\delta$.

# ```math
# \mathbf K \boldsymbol\delta = \mathbf f
# ```

# ## Aeroelasticity

# The vortex lattice method and beam element model are combined into a coupled system to perform static aeroelastic analyses. The analysis is made nonlinear via promotion of the angle of attack $\alpha$ to a variable by specifying the load factor $n$ with a given weight $W$ at fixed sideslip angle $\beta$.

# Define $\mathbf x = [\boldsymbol\Gamma, \boldsymbol\delta, \alpha]$ as the state vector satisfying the residual equations:

# ```math
# \begin{align*}
#     \mathcal R_A(\mathbf x) & = \mathbf A(\boldsymbol\delta) \boldsymbol\Gamma - \mathbf V_\infty(\alpha) \cdot [\mathbf n_i(\boldsymbol\delta)]_{i = 1,\ldots, N} \\
#     \mathcal R_S(\mathbf x) & = \mathbf K \boldsymbol\delta - \mathbf f(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) \\
#     \mathcal R_L(\mathbf x) & = L(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) - n W \\
# \end{align*}
# ```

# where the lift $L$ is obtained from the Kutta-Jowkowski theorem.

# ```math
# [D, Y, L] = \mathbf R_B^W(\alpha, \beta)\left(\sum_{i = 1}^N \rho \mathbf V_i \times  \boldsymbol\Gamma_i \mathbf l_i(\boldsymbol\delta)\right)
# ```

# and the structural load vector $\mathbf f$ is obtained from conservative and consistent load averaging of the nearfield Kutta-Jowkowski forces.

# ```math

# ```

# ## Flight Dynamics

# ### Longitudinal Motion

# A standard three degrees-of-freedom rigid body model is used for performing flight dynamics analyses in the 2-dimensional longitudinal plane formulated as an initial-value problem. The coupled differential equations governing in a canonical state-space representation, with the time-evolution of the state vector $\mathbf x$ driven by the forcing function $\mathbf f$ are shown below subject to initial conditions $\mathbf x_0$.

# ```math
# \begin{align*}
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
# \end{align*}
# ```
# The fuel burn over time (viz. reduction of mass $m$) is computed using a specific fuel consumption value with a linear dependence on the thrust. A manual controller law for the elevator deflection angle $\delta_e$ can also be implemented by providing the function $g(\mathbf x, t)$.

# ### Full Space

# A standard six degrees-of-freedom rigid body model is used for performing flight dynamics analyses in 3 dimensions.

# ```math

# ```

# ## References

# 1. Mark Drela. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
# 2. Joseph Katz and Allen Plotkin. _Low-Speed Aerodynamics, Second Edition_. Cambridge University Press, 2001.
