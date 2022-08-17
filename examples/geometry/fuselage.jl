using AeroFuse

## Fuselage definition
rads = range(0, 1, length = 10)
lens = 0.05:0.05:1.0
fuse = Fuselage(lens[1:end-1], [ rads; reverse(rads) ] )

## Property checks
projected_area(fuse)
length(fuse)

## Cosine interpolation
lens_rads = cosine_interpolation(fuse, 40)

# Circles for plotting
n_pts          = 20
circle3D(r, n) = let arcs = 0:2π/n:2π; [ zeros(length(arcs)) r * cos.(arcs) r * sin.(arcs) ] end

##
xs = lens_rads[:,1]
circs = [ reduce(hcat, eachrow(circ) .+ Ref([x; 0; 0]))' for (x, circ) in zip(xs, circle3D.(lens_rads[:,2], n_pts)) ]

## 
using Plots
pyplot()
plt = plot()
[ plot!(circ[:,1], circ[:,2], circ[:,3], color = :blue, label = :none) for circ in circs ] 
plot!()