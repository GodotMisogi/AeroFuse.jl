# # Theory

# The following theory described is meant to be minimal, with mainly the equations presented for reference when reading the code.

# ## Geometry

# ### Wing Parametrization
# 
# A **half-wing** is defined in terms of a nested trapezoid. A single trapezoid is called a **section**. A section consists of two foil profiles, and their associated chord lengths and twist angles. Between them is their span length with associated _leading-edge_ dihedral and sweep angles. So a general half-wing consisting of ``n`` sections will have ``n`` entries for foils, chords, and twists, and ``n - 1`` entries for spans, dihedrals, sweeps for some ``n \in \mathbb N``. The following illustration should help visualize the concept.
# 
# ![](https://godot-bloggy.xyz/post/diagrams/WingGeometry.svg)

# ## Aerodynamics

# The aerodynamic analyses in AeroFuse mainly utilize potential flow theory and solve problems using a _boundary element method_. This essentially is the following Laplace equation problem with the following Robin (?) boundary conditions:

# ```math
# \nabla^2 \phi = 0, \quad \mathbf V \equiv \nabla \phi \cdot \hat{\mathbf n} = 0, \quad \lim_{\mathbf r \to \infty} \phi(\mathbf r) \to 0
# ```

# !!! note
#     Implementations of viscous-inviscid coupled analyses for drag prediction (á là [XFOIL](https://web.mit.edu/drela/Public/web/xfoil/)) are in progress.

# ### Doublet-Source Panel Method

# The doublet-source panel method predicts **inviscid, incompressible, irrotational, isentropic** external flow over surfaces in 2 dimensions. 

# Source and doublet singularities are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the corresponding free-field Green function for each singularity. 

# ```math

# ```


# The velocities are added to obtain the total induced velocity at a point $\mathbf r$.

# ```math

# ```

# ### Vortex Lattice Method

# The vortex lattice method predicts **inviscid, incompressible, irrotational, isentropic** external flow over "thin" surfaces in 3 dimensions. 

# Vortex filaments are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the Biot-Savart integral for a vortex line of length $\ell$ with a constant circulation strength $\Gamma$.

# ```math
# \mathbf V(\mathbf r, \mathbf r') = \frac{\Gamma}{4\pi} \int_0^\ell \frac{d\boldsymbol\ell' \times (\mathbf r - \mathbf r')}{|\mathbf r - \mathbf r'|^3}
# ```

# The vortices can be set up in various configurations consisting of bound or semi-infinite filaments, commonly in the form of _horseshoes_ or _vortex rings_.

# 1. _Horseshoe elements:_
#    These are defined by a finite _bound leg_ and two semi-infinite _trailing legs_. AeroFuse encodes this information in the `Horseshoe` type.
# 2. _Vortex rings:_
#    These are defined by four _bound legs_. AeroFuse encodes this information in the `VortexRing` type.

# A quasi-steady freestream condition with velocity $\mathbf U$ and rotation $\boldsymbol\Omega$ (in the body's frame) defines an external flow. The induced velocity at a point is given by:

# ```math
# \mathbf V_{\infty}(\mathbf r) = - (\mathbf U + \boldsymbol\Omega \times \mathbf r)
# ```

# The velocities are added to obtain the total induced velocity at a point $\mathbf r$.

# ```math
# \mathbf V(\mathbf r) = \sum_i \frac{\Gamma_i}{4\pi} \int_0^{\ell_i} \frac{d\boldsymbol\ell_i' \times (\mathbf r - \mathbf r_i')}{|\mathbf r - \mathbf r_i'|^3} + \mathbf V_\infty(\mathbf r)
# ```

# Imposing the Neumann boundary condition $\mathbf V \cdot \hat{\mathbf n} = 0$ defines the problem. The construction, in essence, fundamentally results in the following linear system to be solved:

# ```math
# \mathbf A \boldsymbol\Gamma = -V_{\infty} \cdot [\hat{\mathbf n}_i ]_{i = 1, \ldots, N}
# ```

# ### Compressibility Corrections
# The following **Prandtl-Glauert** equation, represented in wind axes, is applied for a weakly compressible flow problem ($0.3 \leq M_\infty \leq 0.7$).

# ```math
# \beta_{PG}^2\frac{\partial^2\phi}{\partial x^2} + \frac{\partial^2\phi}{\partial y^2} + \frac{\partial^2\phi}{\partial z^2} = 0, \quad \beta_{PG}^2 = \left(1 - M_\infty^2\right)
# ```

# The Prandtl-Glauert transformation $\phi(x,y,z; \beta_{PG}) \to \bar\phi(\bar x, \bar y, \bar z)$ converts this equation into an equivalent incompressible flow problem in a transformed geometric space. This "bar" map scales the coordinates $(\bar x,\bar y, \bar z) = (x,\beta_{PG} y, \beta_{PG} z)$ and the potential $\bar\phi = \beta_{PG}^2 \phi$ in sequence. Hence the transformed equation satisfies the Laplace equation with the Neumann boundary condition:
# ```math
# \begin{aligned}
#     \bar\nabla^2 \bar\phi & = 0, \\
#     \bar{\mathbf V} \cdot \hat{\bar{\mathbf n}} & = 0,
# \end{aligned}
# ```
# where $\bar\nabla$ is differentiation with respect to the transformed coordinates and $\hat{\bar{\mathbf n}} = (\beta_{PG} \hat n_x, \hat n_y, \hat n_z)$ which can be proved by computing the appropriate cross product. 

# As the circulation is a scalar (hence invariant of the coordinate transformation but not the potential scaling), the inverse is also readily derived.
# ```math
# \begin{aligned}
#     \Gamma & = \int \mathbf V \cdot d\boldsymbol\ell = \int \nabla\phi \cdot d\boldsymbol\ell, \\
#     \bar{\Gamma} & = \int \bar{\mathbf V} \cdot d\bar {\boldsymbol\ell} = \int\bar\nabla\bar\phi \cdot d\bar{\boldsymbol\ell}, \\
#     \implies \Gamma & = \bar{\Gamma} / \beta_{PG}^2
# \end{aligned}
# ```
# Hence the solution of the resultant incompressible system in transformed coordinates provides the necessary quantities of interest for calculating the dynamics.

# ## Structures

# The structural analyses in AeroFuse utilize _linear finite-element methods_. 

# Particularly, a $2$-dimensional beam element model has been implemented following the standard formulation using cubic Hermite shape functions based on **Euler-Bernoulli beam theory**. These are embedded into a $3$-dimensional local coordinate system in the vortex lattice method without loss of generality.

# The linear system consisting of the stiffness matrix $\mathbf K$ and load vector $\mathbf f$ are solved to obtain the displacement vector $\boldsymbol\delta$.

# ```math
# \mathbf K \boldsymbol\delta = \mathbf f
# ```

# ## Aeroelasticity

# The vortex lattice method and beam element model are combined into a coupled system to perform static aeroelastic analyses. The analysis is made nonlinear via promotion of the angle of attack $\alpha$ to a variable by specifying the load factor $n$ with a given weight $W$ at fixed sideslip angle $\beta$.

# Define $\mathbf x = [\boldsymbol\Gamma, \boldsymbol\delta, \alpha]$ as the state vector satisfying the residual equations:

# ```math
# \begin{aligned}
#     \mathcal R_A(\mathbf x) & = \mathbf A(\boldsymbol\delta) \boldsymbol\Gamma - \mathbf V_\infty(\alpha) \cdot [\mathbf n_i(\boldsymbol\delta)]_{i = 1,\ldots, N} \\
#     \mathcal R_S(\mathbf x) & = \mathbf K \boldsymbol\delta - \mathbf f(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) \\
#     \mathcal R_L(\mathbf x) & = L(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) - n W \\
# \end{aligned}
# ```

# where the lift $L$ is obtained by transforming forces computed via the Kutta-Jowkowski theorem into wind axes.

# ```math
# [D, Y, L] = \mathbf R_B^W(\alpha, \beta)\left(\sum_{i = 1}^N \rho \mathbf V_i \times  \boldsymbol\Gamma_i \boldsymbol \ell_i(\boldsymbol\delta)\right)
# ```

# and the structural load vector $\mathbf f$ is obtained from conservative and consistent load averaging of the Kutta-Jowkowski forces in geometric axes.

# ```math

# ```

# ## References

# 1. Mark Drela. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
# 2. Joseph Katz and Allen Plotkin. _Low-Speed Aerodynamics, Second Edition_. Cambridge University Press, 2001.
