function chop_sections(set1, set2, n :: Integer; spacing = "uniform", direction = 1) 
	if lowercase(spacing) == "cosine"
		space = cosine_spacing(0.5, 1., n + 1)
	elseif lowercase(spacing) == "sine"
		space = sine_spacing(0., 1., (n + 1) * direction)
	else
		space = range(0, stop = 1, length = n + 1)
	end

	[ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

coords_chopper(coords, n, spacing = "cosine", direction = 1) = @views [ chop_sections.(coords[1:end-1], coords[2:end], n; spacing = spacing, direction = direction)...; [coords[end]] ]

chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) ∈ zip(lead, trail) ]

chord_chopper(coords, n) = [ [ weighted_vector(chord[1,:], chord[2,:], μ) for μ ∈ cosine_spacing(0.5, 1., n + 1) ] for chord ∈ coords ]

span_chopper(lead, trail, div, spacing = "cosine", dir = 1) = coords_chopper(lead, div, spacing, dir), coords_chopper(trail, div, spacing, dir)

wing_chopper(lead, trail, span_num, chord_num; span_spacing = "cosine", direction = 1) = chord_chopper(chord_sections(span_chopper(lead, trail, span_num, span_spacing, direction)...), chord_num)

"""
	make_panels(coords)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
function make_panels(coords)
	spanlist = vectarray.(coords)    
	adjacent_sections = zip(spanlist, spanlist[2:end])
	hcat(( Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) ∈ adjacent_sections )...)
end

wetted_area(panels :: Matrix{<: Panel3D}) = sum(panel -> panel_area(panel), panels)