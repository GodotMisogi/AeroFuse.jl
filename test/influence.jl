##
using AeroMDAO
using CoordinateTransformations
using Plots
plotlyjs()

##
panel = Panel3D(
    [1.0, -1.0, 0.0], 
    [0.0, -1.0, 1.0], 
    [0.0, 0.0, 1.0], 
    [1.0, 0.0, 0.0]
)

# Axis permutations
P  = [ 0  1  0 ;
       1  0  0 ;
       0  0 -1 ]

T  = AeroMDAO.PanelGeometry.get_transformation(panel, P)
mv = T.linear[1,:]
lv = T.linear[2,:]
nv = T.linear[3,:]

mp  = midpoint(panel)
mpm = [ mp mp + mv ]' # Streamwise
mpl = [ mp mp + lv ]' # Lateral
mpn = [ mp mp + nv ]' # Normal

##
ε = 1
vec = nv + mv + lv
p = mp + vec * ε
mpp = [ mp inv(T)(T(mp) + T(p)) ]'

@time φ = quadrilateral_doublet_potential(1., panel, p)
@time σ = quadrilateral_source_potential(1., panel, p)

##
pp = plot_panel(panel)
cp = combinedimsview(panel_coordinates(panel), (1))

plot(aspect_ratio = 1, xlim = (-5,5), ylim = (-5,5), zlim = (-5,5), xlabel = "x", ylabel = "y", zlabel = "z")
plot!(pp[:,1], pp[:,2], pp[:,3])

plot!(mpm[:,1], mpm[:,2], mpm[:,3], label = "Stream Vector")
plot!(mpl[:,1], mpl[:,2], mpl[:,3], label = "Lateral Vector")
plot!(mpn[:,1], mpn[:,2], mpn[:,3], label = "Normal Vector")
scatter!(cp[:,1], cp[:,2], cp[:,3], label = "Coordinates")
scatter!([mp.x], [mp.y], [mp.z], label = "Midpoint")
scatter!([p[1]], [p[2]], [p[3]], label = "Evaluation Point")
plot!(mpp[:,1], mpp[:,2], mpp[:,3], label = "Vector")

##
εs = -5:0.05:5
φs = map(εs) do ε
    p = mp + vec * ε
    quadrilateral_doublet_potential(1., panel, p)
    # quadrilateral_source_potential(1., panel, p)
end

vs = combinedimsview([ mp + vec * ε for ε in εs ], (1))

##
plot(εs, φs, xlabel = "z", ylabel = "φ")

##
plot(pp[:,1], pp[:,2], pp[:,3], aspect_ratio = 1, xlim = (-5,5), ylim = (-5,5), zlim = (-5,5),)

scatter!(vs[:,1], vs[:,2], vs[:,3], zcolor = φs, ms = 1)