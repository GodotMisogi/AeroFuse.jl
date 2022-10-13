## Fuselage example
using AeroFuse

# Fuselage parameters
l_fuselage = 18. # Length (m)
w_fuselage = 1.8 # Width (m)

# Chordwise locations and corresponding radii
lens = [0.0, 0.005, 0.01, 0.03, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.98, 1.0]
rads = [0.05, 0.15, 0.25, 0.4, 0.8, 1., 1., 1., 1., 0.85, 0.3, 0.01] * w_fuselage / 2

fuse = Fuselage(l_fuselage, lens, rads, [0., 0., 0.])

## Plotting
using Plots

plot(fuse, aspect_ratio = 1, zlim = (-10,10))

## Aerodynamic model
using StaticArrays
using LinearAlgebra
using Statistics

## Doublet lines
coo = coordinates(fuse, 10)
xs, ys = coo[:,1], coo[:,2]
η = [0.,0.,1.]
coords = SVector.(xs, 0., 0.)
dblt_lines = SourceLine3D.(1.0, coords[1:end-1], coords[2:end])

## Velocity
x = [1.,2.,3.]
vel = velocity.(dblt_lines, Ref(x))

## Panels
fuse_pans = plot_surface(fuse, 10)

# Vectors
ts = diff(coo; dims = 1)
rs = SVector.((xs[1:end-1] + xs[2:end]) / 2, (ys[1:end-1] + ys[2:end]) / 2, 0.)
ns = SVector.(-ts[:,2], ts[:,1], 0.)

## Freestream
fs = Freestream(alpha = 1.0)i
Vinf = velocity(fs) * 5.

## Influence matrix
AIC = [ dot(Vinf + source_velocity(dbl_j, r_i), n_i) for (n_i, r_i) in zip(ns, rs), dbl_j in dblt_lines ]

## Boundary condition
boco = map(n -> -dot(Vinf, n), ns)

##
φs = AIC \ boco

# pans = [ 
#     let coo = fuse_pans[n:(n+1),:,(k+1):-1:k];
#     coords = [ coo[:,:,1]; coo[:,:,2] ];
#     Panel3D(SVector.(coords[:,1], coords[:,2], coords[:,3])...) end;
#     for n in axes(fuse_pans, 1)[1:end-1], k in axes(fuse_pans, 3)[1:end-1] 
# ]