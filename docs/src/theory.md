# Theory

The theory described is meant to be minimal, with mainly the equations presented for reference when reading the code.

## Aerodynamics

The aerodynamic analyses in AeroMDAO mainly utilize potential flow theory and solve problems using a _boundary element method_. This essentially is the following Laplace equation problem with a Neumann boundary condition:

```math
\nabla^2 \phi = 0, \quad \mathbf V \equiv \nabla \phi \cdot \hat{\mathbf n} = 0
```

Implementations for viscous-inviscid analyses are in progress.

### Doublet-Source Panel Method

The doublet-source panel method predicts **inviscid, incompressible, irrotational, isentropic** external flow over surfaces in $2$ dimensions. 

Source and doublet singularities are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the corresponding free-field Green function for each singularity. 

```math

```


As the system is linear, the velocities are added together to obtain the total induced velocity at a point $\mathbf r$.

```math

```

### Vortex Lattice Method

The vortex lattice method predicts **inviscid, incompressible, irrotational, isentropic** external flow over "thin" surfaces in $3$ dimensions. 

Vortex filaments are placed on the surface, and boundary conditions are imposed on their induced velocity to obtain a well-posed problem. The induced velocity is evaluated by the Biot-Savart integral for a vortex line of length $\ell$ with a constant circulation strength $\Gamma$.

```math
\mathbf V(\mathbf r, \mathbf r') = \frac{\Gamma}{4\pi} \int_0^\ell \frac{d\boldsymbol\ell' \times (\mathbf r - \mathbf r')}{|\mathbf r - \mathbf r'|^3}
```

A quasi-steady freestream condition with velocity $\mathbf U$ and rotation $\boldsymbol\Omega$ (in the body's frame) is imposed for the external flow. The induced velocity at a point is given by:

```math
\mathbf V_{\infty}(\mathbf r) = - (\mathbf U + \boldsymbol\Omega \times \mathbf r)
```

As the system is linear, the velocities are added together to obtain the total induced velocity at a point $\mathbf r$.

```math
\mathbf V(\mathbf r) = \sum_i \frac{\Gamma_i}{4\pi} \int_0^{\ell_i} \frac{d\boldsymbol\ell_i' \times (\mathbf r - \mathbf r_i')}{|\mathbf r - \mathbf r_i'|^3} + \mathbf V_\infty(\mathbf r)
```

Imposing the Neumann boundary condition $\mathbf V \cdot \hat{\mathbf n} = 0$ defines the problem. There are various methods for setting up the analysis (using horseshoes or vortex rings), but this is the fundamental construction in essence. results in the following linear system:

```math
\mathbf A \boldsymbol\Gamma = -V_{\infty} \cdot [\hat{\mathbf n}_i ]_{i = 1, \ldots, N}
```

## Structures

The structural analyses in AeroMDAO utilize _linear finite-element methods_. 

Particularly, a $2$-dimensional beam element model has been implemented following the standard formulation using cubic Hermite shape functions based on **Euler-Bernoulli** beam theory. These are embedded into a $3$-dimensional local coordinate system in the vortex lattice method without loss of generality.

The linear system consisting of the stiffness matrix $\mathbf K$ and load vector $\mathbf f$ are solved to obtain the displacement vector $\boldsymbol\delta$.

```math
\mathbf K \boldsymbol\delta = \mathbf f
```

## Aeroelasticity

The vortex lattice method and beam element model are combined into a coupled system to perform static aeroelastic analyses. The analysis is made nonlinear by making the angle of attack $\alpha$ a variable by specifying the load factor $n$ with a given weight $W$ at fixed sideslip angle $\beta$.

Define $\mathbf x = [\boldsymbol\Gamma, \boldsymbol\delta, \alpha]$ as the state vector satisfying the residual equations:

```math
\begin{align*}
    \mathcal R_A(\mathbf x) & = \mathbf A(\boldsymbol\delta) \boldsymbol\Gamma - \mathbf V_\infty(\alpha) \cdot [\mathbf n_i(\boldsymbol\delta)]_{i = 1,\ldots, N} \\
    \mathcal R_S(\mathbf x) & = \mathbf K \boldsymbol\delta - \mathbf f(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) \\
    \mathcal R_L(\mathbf x) & = L(\boldsymbol\delta, \boldsymbol\Gamma, \alpha) - n W \\
\end{align*}
```

where the total lift $L$ is obtained by computing the forces via the Kutta-Jowkowski theorem. The forces are expressed in wind axes by the linear transformation $\mathbf R_G^W$.

```math
[D, Y, L] = \mathbf R_G^W(\alpha, \beta)\left(\sum_{i = 1}^N \rho \mathbf V_i \times \boldsymbol\Gamma_i \boldsymbol \ell_i(\boldsymbol\delta)\right)
```

and the structural load vector $\mathbf f$ is obtained from conservative and consistent load averaging of Kutta-Jowkowski forces:


## Flight Dynamics

A standard six degrees-of-freedom rigid body model is used for performing flight dynamics analyses.

```math

```

## References

1. Mark Drela. _Flight Vehicle Aerodynamics_. MIT Press, 2014.
2. Joseph Katz and Allen Plotkin. _Low-Speed Aerodynamics, Second Edition_. Cambridge University Press, 2001.