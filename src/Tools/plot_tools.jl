plot_panels(panels :: AbstractVector{Panel3D}) = (tupvector ∘ panel_coords).(panels)
plot_streamlines(streams) = tupvector(streams)



## Doublet-source
#==========================================================================================#

## Plotting domain
# x_domain, y_domain = (-1, 2), (-1, 1)
# grid_size = 50
# x_dom, y_dom = linspace(x_domain..., grid_size), linspace(y_domain..., grid_size)
# grid = x_dom × y_dom

# vels, pots = grid_data(dub_src_panels, grid)
# cp = pressure_coefficient.(uniform.mag, vels);

# lower_panels, upper_panels = split_panels(dub_src_panels);

# ## Airfoil plot
# plot( (first ∘ collocation_point).(upper_panels), (last ∘ collocation_point).(upper_panels), 
#         label = "Upper", markershape = :circle,
#         xlabel = "x", ylabel = "C_p")
# plot!((first ∘ collocation_point).(lower_panels), (last ∘ collocation_point).(lower_panels),
#         label = "Lower", markershape = :circle,
#         xlabel = "x", ylabel = "C_p")

# ## Pressure coefficient
# plot( (first ∘ collocation_point).(upper_panels), :cp .<< upper_panels, 
#         label = "Upper", markershape = :circle, 
#         xlabel = "x", ylabel = "C_p")
# plot!((first ∘ collocation_point).(upper_panels), :cp .<< lower_panels, 
#         label = "Lower", markershape = :circle, yaxis = :flip)

# ## Control volume
# p1 = contour(x_dom, y_dom, cp, fill = true)
# plot(p1)
# plot!(first.(:start .<< panels), last.(:start .<< panels), 
#       color = "black", label = "Airfoil", aspect_ratio = :equal)