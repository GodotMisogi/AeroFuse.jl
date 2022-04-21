## Panel helpers
#============================================#

# Source singularities

constant_source_stream(σ, panel :: AbstractPanel2D, xi, yi) = panel_scalar(constant_source_stream, σ, panel, xi, yi)

linear_source_stream_1(σ1, panel :: AbstractPanel2D, xi, yi) =  panel_scalar(linear_source_stream_minus, σ1, panel, xi, yi) - panel_scalar(linear_source_stream_plus, σ1, panel, xi, yi)

linear_source_stream_2(σ2, panel :: AbstractPanel2D, xi, yi) = panel_scalar(linear_source_stream_plus, σ2, panel, xi, yi)

linear_source_stream(σ1, σ2, panel :: AbstractPanel2D, xi, yi) = 
linear_source_stream_1(σ1, panel, xi, yi), linear_source_stream_2(σ2, panel, xi, yi)

source_stream(σ1, σ2, panel :: AbstractPanel2D, xi, yi) = linear_source_stream_1(σ1, panel, xi, yi) + linear_source_stream_2(σ2, panel, xi, yi)

source_stream(σ, panel :: AbstractPanel2D, xi, yi) = panel_scalar(constant_source_stream, σ, panel, xi, yi)

linear_source_velocity(σ1, σ2, panel :: AbstractPanel2D, xi, yi) = panel_velocity(linear_source_velocity_a, σ1, panel, xi, yi), panel_velocity(linear_source_velocity_b, σ2, panel, xi, yi)

source_velocity(σ, panel :: AbstractPanel2D, xi, yi) = panel_velocity(constant_source_velocity, σ, panel, xi, yi)

source_velocity(σ1, σ2, panel :: AbstractPanel2D, xi, yi) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ1, σ2, panel, xi, yi)

source_velocity(σ1, σ2, panel :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(linear_source_velocity_a, linear_source_velocity_b, σ1, σ2, panel, panel_i)

source_velocity(σ, panel :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(constant_source_velocity, σ, panel, panel_i)

# Vortex singularities
linear_vortex_stream_1(γ1, panel :: AbstractPanel2D, xi, yi) =  panel_scalar(linear_vortex_stream_minus, γ1, panel, xi, yi) - panel_scalar(linear_vortex_stream_plus, γ1, panel, xi, yi)

linear_vortex_stream_2(γ2, panel :: AbstractPanel2D, xi, yi) = panel_scalar(linear_vortex_stream_plus, γ2, panel, xi, yi)

linear_vortex_stream(σ1, σ2, panel :: AbstractPanel2D, xi, yi) = 
linear_vortex_stream_1(σ1, panel, xi, yi), linear_vortex_stream_2(σ2, panel, xi, yi)

vortex_stream(γ1, γ2, panel :: AbstractPanel2D, xi, yi) = linear_vortex_stream_1(γ1, panel, xi, yi) + linear_vortex_stream_2(γ2, panel, xi, yi)

vortex_stream(γ, panel :: AbstractPanel2D, xi, yi) = vortex_stream(γ, γ, panel :: AbstractPanel2D, xi, yi) 

vortex_velocity(γ1, γ2, panel :: AbstractPanel2D, x, y) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ1, γ2, panel, x, y)

vortex_velocity(γ, panel :: AbstractPanel2D, x, y) = vortex_velocity(γ, γ, panel :: AbstractPanel2D, x, y)

vortex_velocity(γ1, γ2, panel :: AbstractPanel2D, panel_i :: AbstractPanel2D) = panel_velocity(linear_vortex_velocity_a, linear_vortex_velocity_b, γ1, γ2, panel, panel_i)

linear_vortex_velocity(γ1, γ2, panel :: AbstractPanel2D, xi, yi) = panel_velocity(linear_vortex_velocity_a, γ1, panel, xi, yi), panel_velocity(linear_vortex_velocity_b, γ2, panel, xi, yi)

# Total velocity

total_velocity(σ_1j :: Real, σ_2j :: Real, γ_1j :: Real, γ_2j :: Real, panel_j :: AbstractPanel2D, x, y) = source_velocity(σ_1j, σ_2j, panel_j, x, y) + vortex_velocity(γ_1j, γ_2j, panel_j, x, y)

total_velocity(σ_1j, σ_2j, γ_1j, γ_2j, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) = source_velocity(σ_1j, σ_2j, panel_j, panel_i) .+ vortex_velocity(γ_1j, γ_2j, panel_j, panel_i)

# Influence coefficient computations
neumann_influence_coefficient(velocity_func :: F, angle_func, panel_j :: AbstractPanel2D, panel_i :: AbstractPanel2D) where F <: Function = dot(panel_velocity(velocity_func, 1., panel_j, panel_i), angle_func(panel_i))

# dirichlet_influence_coefficient()

## Matrix construction
#============================================#

function two_point_neumann_matrix(vel_a :: F, vel_b :: G, panels_1, panels_2) where F where G
    N     = min(length(panels_1), length(panels_2))
    inf_a = neumann_influence_matrix(vel_a, normal_vector, panels_1, panels_2)
    inf_b = neumann_influence_matrix(vel_b, normal_vector, panels_1, panels_2)
    [ inf_a zeros(N) ] + [ zeros(N) inf_b ]
end

linear_source_neumann_matrix(panels_1, panels_2) = two_point_neumann_matrix(linear_source_velocity_a, linear_source_velocity_b, panels_1, panels_2)
linear_vortex_neumann_matrix(panels_1, panels_2) = two_point_neumann_matrix(linear_vortex_velocity_a, linear_vortex_velocity_b, panels_1, panels_2)

kutta_condition(panels_1, panels_2) = [ 1 zeros(length(panels_1) - 1)' 1 zeros(length(panels_2) - length(panels_1))' ]
source_influence_matrix(panels_1, panels_2) = linear_source_matrix(panels_1, panels_2)
vortex_influence_matrix(panels_1, panels_2) = linear_vortex_matrix(panels_1, panels_2)

# Diagonal cases
kutta_condition(panels) = [ 1 zeros(length(panels) - 1)' 1 ]
source_influence_matrix(panels) = [ linear_source_matrix(panels, panels); kutta_condition(panels) ]
vortex_influence_matrix(panels) = [ linear_vortex_matrix(panels, panels); kutta_condition(panels) ]


constant_source_influence_coefficient(panel_j, panel_i) = ifelse(panel_i === panel_j, 0.5, influence_coefficient(constant_source_velocity, normal_vector, panel_j, panel_i))

constant_source_matrix(panels) = [ constant_source_influence_coefficient(panel_j, panel_i) for (panel_j, panel_i) in product(panels, panels)]

neumann_influence_matrix(func, angle_func, panel_is, panel_js) = [ influence_coefficient(func, angle_func, panel_j, panel_i) for (panel_i, panel_j) in product(panel_is, panel_js) ]

function two_point_dirichlet_matrix(src_a :: F, src_b :: G, panels_1, panels_2) where {F,G}
    N     = min(length(panels_1), length(panels_2))
    inf_a = dirichlet_influence_matrix(src_a, normal_vector, panels_1, panels_2)
    inf_b = dirichlet_influence_matrix(src_b, normal_vector, panels_1, panels_2)
    [ inf_a zeros(N) ] + [ zeros(N) inf_b ]
end


# Boundary conditions
#============================================#

constant_source_boundary_condition(panels, u) = -map(pan -> dot(u, normal_vector(pan)), panels)

neumann_boundary_condition(panels, u) = [ constant_source_boundary_condition(panels, u); 0 ]


## Block matrix setup
#===============================================#

# linear_vortex_matrix()

function influence_matrix(panels, pts)
    N   = length(pts[:,1])

    # Linear vorticity-streamfunction matrices
    A1 = [ linear_vortex_stream_1(1., pj, xi[1], xi[2]) for xi in pts, pj in panels ]
    A2 = [ linear_vortex_stream_2(1., pj, xi[1], xi[2]) for xi in pts, pj in panels ]
    
    AIC_pts = [ A1 zeros(N) ] + [ zeros(N) A2 ]

    ## Trailing edge
    te_panel, h_TE, tcp, tdp, _ = trailing_edge_info(panels)

    # Constant-strength vortex and source terms
    A_TE = [ 
            linear_vortex_stream_1(1., te_panel, xi[1], xi[2]) + linear_vortex_stream_2(1., te_panel, xi[1], xi[2]) 
            for xi in pts
            ]
    B_TE = [
            constant_source_stream(1., te_panel, xi[1], xi[2]) for xi in pts ] 

    AB_TE = (A_TE * tdp - B_TE * tcp) / 2

    AIC_TE = [ AB_TE zeros(N, N - 2) -AB_TE ] 

    AIC_surf = AIC_pts + AIC_TE

    # Check trailing edge thickness and augment last equation
    if abs(h_TE) < 1e-10 
        AIC[end-1,:] .= 0
        AIC[end-1,[1,2,3,end-3,end-2,end-1]] = [1,-2,1,-1,2,-1]
    end

    # Kutta condition
    kutta = [ 1. zeros(1, N - 2) 1. 0. ]

    # Assembly
    AIC  = [ AIC_surf  fill(-1., N);
               kutta            ]

    AIC
end

boundary_condition(pts) = @views [ -getindex.(pts, 2) getindex.(pts, 1);
                                    zeros(1, 2) ]


function solve_linear(panels)
    # Panel points
    pts  = panel_points(panels)
    AIC  = influence_matrix(panels, pts)
    boco = boundary_condition(pts)
    ψs   = AIC \ boco

    ψs, AIC, boco
end

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
    te_panel, h_TE, tcp, tdp, _ = trailing_edge_info(panels)

    # Construct everything in one loop
    for i in eachindex(pts)
        x, y = pts[i]
        for j in eachindex(panels)
            AIC[i,j]   += linear_vortex_stream_1(1., panels[j], x, y)
            AIC[i,j+1] += linear_vortex_stream_2(1., panels[j], x, y)
        end
        
        # Trailing edge influence coefficients

        # Constant-strength vortex panel
        a_TE = linear_vortex_stream_1(1., te_panel, x, y) + linear_vortex_stream_2(1., te_panel, x, y)
        AIC[i,1]     += a_TE * ( tdp / 2)
        AIC[i,end-1] += a_TE * (-tdp / 2)

        # Constant-strength source panel
        b_TE = constant_source_stream(1., te_panel, x, y)
        AIC[i,1]     += b_TE * (-tcp / 2)
        AIC[i,end-1] += b_TE * ( tcp / 2)

        # Constant streamfunction value on the surface
        AIC[i,end] = -1

        # Check trailing edge thickness and augment last equation accordingly
        if abs(h_TE) < 1e-10 
            AIC[end-1,:] .= 0
            AIC[end-1,[1,2,3,end-3,end-2,end-1]] = [1,-2,1,-1,2,-1]
        end

        # RHS
        boco[i,1] = -y
        boco[i,2] =  x
    end
    
    # Kutta condition
    AIC[end,1]     = 1
    AIC[end,end-1] = 1

    nothing
end
