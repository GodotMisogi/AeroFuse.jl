## Example script for testing linear-strength vortex and source panel methods
using LinearAlgebra
using Base.Iterators
using Seaborn
using AeroMDAO

## Tests
wedge_pts 	 = @. Point2D([0.5, 0.0, -0.5, 0.0, 0.5],
                          [-0.0, -0.5, 0.0, 0.5, 0.0])
wedge_panels = @. Panel2D(wedge_pts[1:end-1], wedge_pts[2:end])
uniform 	 = Uniform2D(1.0, 0.0)

## Source
A = source_influence_matrix(wedge_panels)
b = neumann_boundary_condition(wedge_panels, velocity(uniform))
σs = A \ b

##
sum(σs)

##
x_domain, y_domain = (-1, 1), (-1, 1)
grid_size = 50
x_dom, y_dom = range(x_domain..., length = grid_size), range(y_domain..., length = grid_size)
grid = product(x_dom, y_dom)
wedge = panel_points(wedge_panels)

##
# 
source_vels = [ velocity(uniform) .+ sum(source_velocity.(σs[1:end-1], σs[2:end], wedge_panels, x, y)) for (x, y) in grid ]
# source_vels = [ sum(source_velocity.(1., 1., wedge_panels, x, y)) for (x, y) in grid ]
speeds 		= @. sqrt(first(source_vels)^2 + last(source_vels)^2)

cock = plot(first.(wedge), last.(wedge))
streamplot(first.(grid)', last.(grid)', first.(source_vels)', last.(source_vels)', color = speeds', density = 2.5)

show()


## Vortex
A = vortex_influence_matrix(wedge_panels)
b = neumann_boundary_condition(wedge_panels, velocity(uniform))
γs = A \ b

##
sum(γs)

##
vortex_vels = [ sum(vortex_velocity(1., 1., panel, x, y) for panel in wedge_panels) for (x, y) in grid ]

vortex_vels = [ velocity(uniform) .+ sum(vortex_velocity.(γs[1:end-1], γs[2:end], wedge_panels, x, y)) for (x, y) in grid ]

speeds 		= @. sqrt(first(vortex_vels)^2 .+ last(vortex_vels)^2)

plot(first.(wedge), last.(wedge))
streamplot(first.(grid)', last.(grid)', first.(vortex_vels)', last.(vortex_vels)', color = speeds', density = 2.5)

show()

##
vortex_source_vels = [ sum(total_velocity.(σs[1:end-1], σs[2:end], γs[1:end-1], γs[2:end], wedge_panels, x, y)) for (x, y) in grid ]
speeds 			   = @. sqrt(first(vortex_source_vels)^2 + last(vortex_source_vels)^2)

plot(first.(wedge), last.(wedge))
quiver(first.(grid), last.(grid), first.(vortex_source_vels), last.(vortex_source_vels), speeds)

show()

