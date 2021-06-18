function chop_sections(set1, set2, n :: Integer, spacing = "cosine"; flip = false)
	if lowercase(spacing) == "uniform"
		space = uniform_spacing(0., 1., n + 1)
	elseif lowercase(spacing) == "sine"
		space = sine_spacing(0., 1., (n + 1) * ifelse(flip, -1, 1))
	else
		space = cosine_spacing(0.5, 1., n + 1)
	end

	[ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

chop_coordinates(coords, n, spacings = "cosine", flip = false) = @views [ chop_sections.(coords[1:end-1], coords[2:end], n, spacings; flip = flip)...; [coords[end]] ]

chord_sections(lead, trail) = [ [ l'; t' ] for (l, t) ∈ zip(lead, trail) ]

chop_chords(coords, n) = @views [ [ weighted_vector(chord[1,:], chord[2,:], μ) for μ ∈ cosine_spacing(0.5, 1., n + 1) ] for chord ∈ coords ]

chop_spans(lead, trail, div, spacing = "cosine", flip = false) = chop_coordinates(lead, div, spacing, flip), chop_coordinates(trail, div, spacing, flip)

chop_wing(lead, trail, span_num, chord_num; span_spacing = "cosine", flip = false) = chop_chords(chord_sections(chop_spans(lead, trail, span_num, span_spacing, flip)...), chord_num)

"""
	make_panels(xyzs)

Convert an array of coordinates corresponding to a wing, ordered from root to tip and leading-edge to trailing-edge, into panels.
"""
make_panels(xyzs) = @views Panel3D.(xyzs[1:end-1,1:end-1], xyzs[1:end-1,2:end], xyzs[2:end,1:end-1], xyzs[2:end,2:end])
	
# WTF was I thinking?
# spanlist = vectarray.(coords)
# spanlist = zip(coords, coords[2:end,:])
# adjacent_sections = zip(spanlist, spanlist[2:end])
# @views hcat(( Panel3D.(root[1:end-1], root[2:end], tip[2:end], tip[1:end-1]) for (root, tip) ∈ adjacent_sections )...)