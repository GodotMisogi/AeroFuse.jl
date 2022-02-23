## Matrix helpers
#===========================================================================#

linear_vortex_stream_1(γ1, panel_i :: AbstractPanel2D, point) =  panel_scalar(linear_vortex_stream_minus, γ1, panel_i, first(point), last(point)) - panel_scalar(linear_vortex_stream_plus, γ1, panel_i, first(point), last(point))

linear_vortex_stream_2(γ2, panel_i :: AbstractPanel2D, point) = panel_scalar(linear_vortex_stream_plus, γ2, panel_i, first(point), last(point))

constant_source_stream(σ, panel_i :: AbstractPanel2D, point) = panel_scalar(constant_source_stream, σ, panel_i, first(point), last(point))

## Preallocated setup
#===============================================#

function solve_linear!(A, b, panels :: AbstractVector{<: Panel2D{T}}) where T <: Real
    assemble_system!(A, b, panels)
    ψ = A \ b
end

function assemble_system!(AIC, boco, panels)
    # Panel points
    pts = panel_points(panels)

    # Trailing edge info
    te_panel  = trailing_edge_panel(panels)
    s    = panel_vector(reverse_panel(te_panel))
    t    = normalize(bisector(panels[1], panels[end]))
    p    = normalize(s)
    h_TE = det([ permutedims(t)
                 permutedims(s) ])
    tcp  = norm(t × p)
    tdp  = p ⋅ t

    for i in eachindex(pts)
        xi = pts[i]
        for j in eachindex(panels)
            AIC[i,j]   += linear_vortex_stream_1(1., panels[j], xi)
            AIC[i,j+1] += linear_vortex_stream_2(1., panels[j], xi)
        end
        
        # Trailing edge influence coefficients

        # Constant-strength source panel
        a = constant_source_stream(1., te_panel, xi)
        AIC[i,1]     -= a * (tcp / 2)
        AIC[i,end-1] += a * (tcp / 2)

        # Constant-strength vortex panel
        a = linear_vortex_stream_1(1., te_panel, xi) + linear_vortex_stream_2(1., te_panel, xi)
        AIC[i,1]     -= a * (-tdp / 2)
        AIC[i,end-1] += a * (-tdp / 2)

        # Constant streamfunction value on the surface
        AIC[i,end] = -1

        # Check trailing edge thickness and augment last equation
        if abs(h_TE) < 1e-10 
            AIC[end-1,:] .= 0
            AIC[end-1,[1,2,3,end-3,end-2,end-1]] = [1,-2,1,-1,2,-1]
        end

        # RHS (correct)
        boco[i,1] = -xi[2]
        boco[i,2] =  xi[1]
    end
    
    # Kutta condition (correct)
    AIC[end,1] = 1
    AIC[end,end-1] = 1

    nothing
end

# function aerodynamic_coefficients(panels :: Vector{<: AbstractPanel2D}, γs, u)
#       speed = norm(u)
    
#       vels  = sum(γs)
#       Δrs   = midpair_map(distance, panels)
#       cps   = @. pressure_coefficient(speed, vels)
#       cls   = @. lift_coefficient(cps, Δrs / 2, panel_angle(panels))
    
#       cps, cls
# end

# source_matrix(panels_1 :: Vector{<: AbstractPanel2D}, panels_2 :: Vector{<: AbstractPanel2D}) = [ source_influence(panel_j, panel_i) for panel_i ∈ panels_1, panel_j ∈ panels_2 ]

# source_strengths(panels :: Vector{<: AbstractPanel2D}, u) = dot.(Ref(u), normal_vector.(panels))

# boundary_vector(panels :: Vector{<: AbstractPanel2D}, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

# function boundary_vector(panels :: Vector{<: AbstractPanel2D}, wakes :: Vector{<: AbstractPanel2D}, u) 
#       source_panels = [ panels; wakes ]
#       [ - source_matrix(panels, source_panels) * source_strengths(source_panels, u); 0 ]
# end


## Matrix construction
#===========================================================================#

kutta_condition(points) = [ 1 zeros(length(points) - 2)' 1 ]

bisector(f_pan, l_pan) = (normalize(panel_vector(l_pan)) + normalize(-panel_vector(f_pan))) / 2

trailing_edge_influence(s, t, te_panel, pt) = ifelse(p1(te_panel) == pt, 0.0, trailing_edge_source(1., te_panel, pt) * norm(t × s) + ψp(1., 1., te_panel, pt) * (abs ∘ dot)(t, s)) 

function trailing_edge(panels :: Vector{<: AbstractPanel2D}, pts)
    te_panel  = trailing_edge_panel(panels)
    s         = p2(te_panel) - p1(te_panel)
    t         = bisector(panels[1], panels[end])
    te_vector = trailing_edge_influence.(Ref(s), Ref(t), Ref(te_panel), pts)

    te_panel, te_vector
end

function influence_matrix(panels :: Vector{<: AbstractPanel2D})
    rev_panels  = reverse_panel.(panels)[end:-1:1]
    pts         = panel_points(rev_panels)
    num_pts     = length(pts)

    ψm = [ ψm_influence_matrix(rev_panels, pts) zeros(num_pts) ] 
    ψp = [ zeros(num_pts) ψp_influence_matrix(rev_panels, pts) ]

    ψ_matrix = ψm + ψp

    # Trailing edge setup
    te_panel, te_vector = trailing_edge(rev_panels, pts)
    te_matrix           = [ te_vector zeros(num_pts, num_pts - 2) -te_vector ]
    r_te                = panel_length(te_panel)
    
    inf_matrix          = ψ_matrix + te_matrix

    ψ0_vector           = fill(-1, length(pts))

    [     inf_matrix        ψ0_vector  ;
      kutta_condition(pts)      0      ]
end

boundary(colpoints, uni) = - [ [ stream(uni, x, y) for (x, y) in colpoints ]; 0 ]
