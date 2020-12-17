using PlotlyJS

plot_panels(panels :: AbstractVector{Panel3D}) = panel_coords.(panels)

panel_splits(coords) = [ [ c[1] for c in panel ] for panel in coords ], [ [ c[2] for c in panel ] for panel in coords ], [ [ c[3] for c in panel ] for panel in coords ]

trace_coords(xs, ys, zs, color = :black) = [ scatter3d(
                                                x = x,
                                                y = y,
                                                z = z,
                                                mode = :lines, 
                                                line = attr(color = color),
                                                showlegend = false,
                                                ) 
                                            for (x, y, z) in zip(xs, ys, zs) ]

function trace_panels(panels :: AbstractVector{Panel3D})
    coords = plot_panels(panels)
    xs, ys, zs = panel_splits(coords)
    trace_coords(xs, ys, zs)
end

trace_surface(wing :: Union{HalfWing, Wing}, span_num = 5, chord_num = 30; rotation = one(RotMatrix{3, Float64}), translation = SVector(0, 0, 0)) = (trace_panels ∘ collect)(transform(panel, rotation, translation) for panel in mesh_wing(wing, span_num, chord_num)[:])

function trace_panels(panels :: AbstractVector{Panel3D}, Γs :: AbstractVector{<: Real})
    coords = plot_panels(panels[:])

    min_Γ, max_Γ = extrema(Γs)
    Γ_range = -map(-, min_Γ, max_Γ)
    norm_Γs = [ 2 * (Γ - min_Γ) / Γ_range - 1 for Γ ∈ Γs ];

    xs, ys, zs = panel_splits(coords)
    
    trace = [ mesh3d(
                     x = x,
                     y = y,
                     z = z,
                     intensity = fill(norm_Γ, length(x)),
                     text = norm_Γ,
                     showscale = false,
                     ) 
                for (x, y, z, norm_Γ) in zip(xs, ys, zs, norm_Γs) ]
end

function trace_streamlines(freestream :: Freestream, horseshoe_panels :: AbstractVector{Panel3D}, horseshoes :: AbstractVector{Horseshoe}, Γs :: AbstractVector{<: Real}, length :: Real, num_steps :: Integer = 100)
    streams = streamlines(freestream, horseshoe_panels[:], horseshoes, Γs, length, num_steps)
    xs, ys, zs = panel_splits(streams)
    trace_coords(xs, ys, zs, :lightblue)
end


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