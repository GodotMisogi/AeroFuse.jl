plot_panels(panels :: Array{Panel3D}) = (tupvector ∘ panel_coords).(panels)
plot_streamlines(streams) = tupvector(streams)