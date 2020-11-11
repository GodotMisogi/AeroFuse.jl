##
# coords_xs = [ [ c[1] for c in panel ] for panel in horseshoe_coords ]
# coords_ys = [ [ c[2] for c in panel ] for panel in horseshoe_coords ]
# coords_zs = [ [ c[3] for c in panel ] for panel in horseshoe_coords ]
# streams_xs = [ [ c[1] for c in panel ] for panel in streams ]
# streams_ys = [ [ c[2] for c in panel ] for panel in streams ]
# streams_zs = [ [ c[3] for c in panel ] for panel in streams ]


##
# layout = Layout(;
#                 title = "Penguins",
#                 xaxis_title = "x",
#                 yaxis_title = "y",
#                 aspectratio = attr(x=1, y=15, z=15)
#                 )

# trace_coords = [ mesh3d(
#                         x = x,
#                         y = y,
#                         z = z,
#                         # colorscale = "rainbow",
#                         intensity = repeat([norm_Γ], length(x)),
#                         showscale = false,
#                         ) for (x, y, z, norm_Γ) in zip(coords_xs, coords_ys, coords_zs, norm_Γs) ]

# trace_streams = [ scatter3d(
#                             x = x, 
#                             y = y, 
#                             z = z, 
#                             mode = :lines, 
#                             line = attr(color =:green),
#                             showlegend = false,
#                             ) for (x, y, z) in zip(streams_xs, streams_ys, streams_zs) ]

# plot([ 
#         [ trace for trace in trace_coords ]..., 
#         [ trace for trace in trace_streams ]... 
#      ], 
#      layout)