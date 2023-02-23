## Fuselage example
using AeroFuse

# Fuselage parameters
l_fuselage = 18. # Length (m)
w_fuselage = 1.8 # Width (m)

# Chordwise locations and corresponding radii
lens = [0.0, 0.005, 0.01, 0.03, 0.1, 0.2, 0.4, 0.6, 0.7, 0.8, 0.98, 1.0]
rads = [0.05, 0.15, 0.25, 0.4, 0.8, 1., 1., 1., 1., 0.85, 0.3, 0.01] * w_fuselage / 2

fuse = Fuselage(l_fuselage, lens, rads, [0., 0., 0.])

## Aerodynamic model
using StaticArrays
using LinearAlgebra
using Statistics
using StructArrays

## Doublet lines
n = 11 # Number of points
coo = coordinates(fuse, n)
xs, ys = coo[:,1], coo[:,2]
coords = SVector.(xs, 0., 0.)
src_lines = SourceLine3D.(1.0, coords[1:end-1], coords[2:end])
dblt_lines = DoubletLine3D.(1.0, coords[1:end-1], coords[2:end], Ref(@SVector [0., 0., 1.]))

## Velocity
x = [1.,2.,3.]
vel = velocity.(src_lines, Ref(x))

## Panels
fuse_pans = plot_surface(fuse, 10)

# Vectors
ts = diff(coo; dims = 1)
rs = SVector.((xs[1:end-1] + xs[2:end]) / 2, (ys[1:end-1] + ys[2:end]) / 2, 0.)
ns = SVector.(-ts[:,2], ts[:,1], 0.)

## Freestream
fs = Freestream(alpha = 1.0)
Vinf = velocity(fs) * 5.

## Influence matrix
AIC = [ 
    # potential(fs, r_i) + potential(line_j, r_i) # Dirichlet
    dot( # Neumann
        Vinf + velocity(line_j, r_i),
        n_i
    ) 
    for (n_i, r_i) in zip(ns, rs), line_j in src_lines 
]

## Boundary condition
boco = map(n -> -dot(Vinf, n), ns)

##
σs = AIC \ boco

using Accessors
src_lines = StructArray(map((line, σ) -> line(strength = σ), src_lines, σs))



## Streamlines
seed = cfs.(Spherical(0.001, 0.0, θ) for θ in LinRange(0, 2π, 51)) .- Ref([-0., 0., 0.])
tran = CoordinateTransformations.LinearMap(AngleAxis(π/2, 1, 0., -1))
seed = tran.(seed)

## Grid
Nt = 50
Np = length(seed)
Δt = 0.1
L  = length(fuse)

# Streamlines
Ps = zeros(Nt, Np, 3)

# Initial condition point
Ps[1,:,:] = combinedimsview(seed)'

# Iterate over seed
for (i, _) in enumerate(seed)
    # Iterate over length-steps
    for j in 1:Nt-1
        # Compute new velocity
        V = Vinf + sum(line -> velocity(line, Ps[j,i,:]), src_lines)

        # Compute next point
        Ps[j+1,i,:] = Ps[j,i,:] + V / norm(V) * L / Np
    end
end


##
ls = Point3f.(splitdimsview(Ps, (1,2)))


## Visualization
using WGLMakie
WGLMakie.activate!()

fig = Figure()
scene = LScene(fig[1,1])

[ lines!(stream, color = :green) for stream in eachcol(ls) ]

lines!(Point3f.(rs))
arrows!(Point3f.(rs), Vec3f.(ns), 
    arrowsize = Vec3f(0.3, 0.3, 0.4),
    fxaa = true,
    normalize = true,
    align = :center
)

fig
## Plotting
# using Plots

# plot(fuse, aspect_ratio = 1, zlim = (-10,10))


# pans = [ 
#     let coo = fuse_pans[n:(n+1),:,(k+1):-1:k];
#     coords = [ coo[:,:,1]; coo[:,:,2] ];
#     Panel3D(SVector.(coords[:,1], coords[:,2], coords[:,3])...) end;
#     for n in axes(fuse_pans, 1)[1:end-1], k in axes(fuse_pans, 3)[1:end-1] 
# ]

