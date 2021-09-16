## Matrix helpers
#===========================================================================#

function wake_panels(panels, bound, num)
	lastx, lasty = (p2 ∘ last)(panels)
	bounds = range(lastx, bound, length = num)
	@. WakePanel2D(SVector(bounds[1:end-1], lasty), SVector(bounds[2:end], lasty))
end

ψp(str1, str2, panel_i :: AbstractPanel2D, point) = panel_scalar(vortex_stream_plus, str1, panel_i, first(point), last(point)) + panel_scalar(vortex_stream_minus, str2, panel_i, first(point), last(point))

ψm(str1, str2, panel_i :: AbstractPanel2D, point) = panel_scalar(vortex_stream_plus, str1, panel_i, first(point), last(point)) - panel_scalar(vortex_stream_minus, str2, panel_i, first(point), last(point))

trailing_edge(str1, str2, panel_i :: AbstractPanel2D, point) = panel_scalar(source_stream, str1, panel_i, first(point), last(point)) + panel_scalar(vortex_stream_plus, str2, panel_i, first(point), last(point))

source(str, panel_i, point) = panel_scalar(source_stream, str, panel_i, first(point), last(point))

trailing_edge_source(str, panel_i, point) = panel_scalar(source_stream_trailing_edge, str, panel_i, first(point), last(point))

## Matrix construction
#===========================================================================#

# Need to check singularity conditions
ψp_matrix(panels :: Vector{<: AbstractPanel2D}, points) = [ ψp(1., 1., panel_j, point_i) for (point_i, panel_j) in product(points, panels) ]

ψm_matrix(panels :: Vector{<: AbstractPanel2D}, points) = [ ψm(1., 1., panel_j, point_i) for (point_i, panel_j) in product(points, panels) ]

kutta_condition(points) = [ 1 zeros(length(points) - 2)' 1 ]

streamline(f_pan, l_pan) = (p2(l_pan) - p1(l_pan) + p1(f_pan) - p2(f_pan)) / 2

trailing_edge_influence(s, t, te_panel, pt) = ifelse(p1(te_panel) == pt, 0.0, trailing_edge_source(1., te_panel, pt) * norm(t × s) + ψp(1., 1., te_panel, pt) * (abs ∘ dot)(t, s)) 

function trailing_edge(panels :: Vector{<: AbstractPanel2D}, pts)
	te_panel  = trailing_edge_panel(panels)
	s 		  = p2(te_panel) - p1(te_panel)
	t 	  	  = streamline(panels[1], panels[end])
	te_vector = trailing_edge_influence.(Ref(s), Ref(t), Ref(te_panel), pts)

	te_panel, te_vector
end

function influence_matrix(panels :: Vector{<: AbstractPanel2D})
	rev_panels 	= reverse_panel.(panels)[end:-1:1]
	pts 		= panel_points(rev_panels)
	num_pts 	= length(pts)

	ψM 			= [ ψm_matrix(rev_panels, pts) zeros(num_pts) ] 
	ψP 			= [ zeros(num_pts) ψp_matrix(rev_panels, pts) ]
	ψM, ψP

	ψ_matrix 	= ψM .+ ψP

	# Trailing edge setup
	te_panel, te_vector = trailing_edge(rev_panels, pts)
	te_matrix 	= [ te_vector zeros(num_pts, num_pts - 2) -te_vector ]
	r_te 		= panel_length(te_panel)
	
	inf_matrix 	= ψ_matrix + te_matrix

	ψ0_vector 	= fill(-1, length(pts))

	[     inf_matrix		ψ0_vector ;
	  kutta_condition(pts)		0 	  ]
end

boundary(colpoints, uni) = - [ [ stream(uni, x, y) for (x, y) in colpoints ]; 0 ]


# function aerodynamic_coefficients(panels :: Vector{<: AbstractPanel2D}, γs, u)
# 	speed = norm(u)
	
# 	vels  = sum(γs)
# 	Δrs   = midpair_map(panel_dist, panels)
# 	cps   = @. pressure_coefficient(speed, vels)
# 	cls   = @. lift_coefficient(cps, Δrs / 2, panel_angle(panels))
	
# 	cps, cls
# end

# source_matrix(panels_1 :: Vector{<: AbstractPanel2D}, panels_2 :: Vector{<: AbstractPanel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

# source_strengths(panels :: Vector{<: AbstractPanel2D}, u) = dot.(Ref(u), panel_normal.(panels))

# boundary_vector(panels :: Vector{<: AbstractPanel2D}, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

# function boundary_vector(panels :: Vector{<: AbstractPanel2D}, wakes :: Vector{<: AbstractPanel2D}, u) 
# 	source_panels = [ panels; wakes ]
# 	[ - source_matrix(panels, source_panels) * source_strengths(source_panels, u); 0 ]
# end
