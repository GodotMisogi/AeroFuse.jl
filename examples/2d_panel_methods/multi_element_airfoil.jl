## Multi-element airfoil using constant-strength doublet-source panel method
using AeroMDAO
using CoordinateTransformations, Rotations

## Geometry variables
alpha_u   = [0.2, 0.3, 0.2, 0.15, 0.2]
alpha_l   = [-0.2, -0.1, -0.1, -0.001]
dzs       = (0., 0.)

## Airfoils
airfoil_1 = kulfan_CST(alpha_u, alpha_l, dzs, (0., 0.), 60);      # Kulfan CST airfoil
airfoil_2 = affine(scale(naca4((0,0,1,2)), 0.2); vector = [1.0, -0.1], angle = deg2rad(5.))  # NACA 4-digit airfoil

## Case setup
panels  = ComponentVector(
                            kulfan = make_panels(airfoil_1, 60),
                            naca   = make_panels(airfoil_2, 60)
                          )
uniform = Uniform2D(1., 0.)



## Obviously wrong because Kutta condition is hand-coded -_-
system  = @time solve_system(
                     panels, uniform, 15, 1e5
                    );

##
@time ues = surface_velocities(system);
@time cl  = lift_coefficient(system)
@time cls, cms, cps = surface_coefficients(system)

##
using CairoMakie

af_1 = coordinates(airfoil_1)
af_2 = coordinates(airfoil_2)

f = Figure()
Axis(f[1,1])
lines!(af_1[:,1], af_1[:,2])
lines!(af_2[:,1], af_2[:,2])
f