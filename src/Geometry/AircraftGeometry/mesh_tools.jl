function chop_sections(set1, set2, n :: Integer; spacing = "uniform") 
	if spacing == "cosine"
		space = cosine_dist(0.5, 1., n + 1)
	else
		space = range(0, stop = 1, length = n + 1)
	end

	[ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

coords_chopper(coords, n, spacing = "cosine") = [ chop_sections.(coords[1:end-1], coords[2:end], n; spacing = spacing)...; [coords[end]] ]

chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) ∈ zip(lead, trail) ]

chord_chopper(coords, n) = [ [ weighted_vector(chord[1,:], chord[2,:], μ) for μ ∈ cosine_dist(0.5, 1., n + 1) ] for chord ∈ coords ]

span_chopper(lead, trail, div, spacing = "cosine") = coords_chopper(lead, div, spacing), coords_chopper(trail, div, spacing)

wing_chopper(lead, trail, span_num, chord_num) = chord_chopper(chord_sections(span_chopper(lead, trail, span_num)...), chord_num)

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