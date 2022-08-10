"""
    doublet_matrix(panels_1, panels_2)

Create the matrix of doublet potential influence coefficients between pairs of `panels₁` and `panels₂`.
"""
doublet_matrix(panels_1, panels_2) = [ doublet_influence(panel_j, panel_i) for panel_i in panels_1, panel_j in panels_2 ]

"""
    doublet_matrix(panels_1, panels_2)

Create the matrix of source potential influence coefficients between pairs of `panels₁` and `panels₂`.
"""
source_matrix(panels_1, panels_2) = [ source_influence(panel_j, panel_i) for panel_i in panels_1, panel_j in panels_2 ]

"""
    kutta_condition(panels)

Create the vector describing Morino's Kutta condition given `Panel2Ds`.
"""
kutta_condition(panels :: AbstractVector{<:AbstractPanel2D}) = [ 1 zeros(length(panels) - 2)' -1 ]

kutta_condition(Nf, Nw) = Matrix([ I(Nw) zeros(Nw, Nf-2*Nw) -I(Nw) -I(Nw) ])

"""
    wake_vector(woke_panel :: AbstractPanel2D, panels)

Create the vector of doublet potential influence coefficients from the wake on the panels given the wake panel and the array of `Panel2Ds`.
"""
wake_vector(woke_panel :: AbstractPanel2D, panels) = doublet_influence.(Ref(woke_panel), panels)

"""
    influence_matrix(panels, wake_panel :: AbstractPanel2D)

Assemble the Aerodynamic Influence Coefficient matrix consisting of the doublet matrix, wake vector, Kutta condition given `Panel2Ds` and the wake panel.
"""
influence_matrix(panels, woke_panel :: AbstractPanel2D) =
    [ doublet_matrix(panels, panels)  wake_vector(woke_panel, panels) ;
        kutta_condition(panels)                      1.               ]

"""
    source_strengths(panels, freestream)

Create the vector of source strengths for the Dirichlet boundary condition ``σ = \\vec U_{\\infty} ⋅ n̂`` given Panel2Ds and a Uniform2D.
"""
source_strengths(panels, u) = dot.(Ref(u), normal_vector.(panels))

"""
    boundary_vector(panels, u)

Create the vector for the boundary condition of the problem given an array of Panel2Ds and velocity ``u``.
"""
boundary_vector(panels, u) = [ - source_matrix(panels, panels) * source_strengths(panels, u); 0 ]

# boundary_vector(colpoints, u, r_te) = [ dot.(colpoints, Ref(u)); dot(u, r_te) ]

boundary_vector(panels, u, r_te) = [ dot.(collocation_point.(panels), Ref(u)); dot(u, r_te) ]

function boundary_vector(panels :: Vector{<: AbstractPanel2D}, wakes :: Vector{<: AbstractPanel2D}, u)
    source_panels = [ panels; wakes ]
    [ -source_matrix(panels, source_panels) * source_strengths(source_panels, u); 0 ]
end

function boundary_vector(panels :: AbstractArray{<: AbstractPanel3D}, wakes, V∞)
    panelview = @view permutedims(panels)[:]
    B = source_matrix(panelview, panelview)
    σ = source_strengths(panelview, V∞)

    return [ -B * σ; zeros(length(wakes)) ]
end

"""
    solve_linear(panels, u, sources, bound)

Solve the system of equations ``[AIC][φ] = [\\vec{U} ⋅ n̂] - B[σ]`` condition given the array of Panel2Ds, a velocity ``\\vec U``, a condition whether to disable source terms (``σ = 0``), and an optional named bound for the length of the wake.
"""
function solve_linear(panels, u, α, r_te, sources :: Bool; bound = 1e2)
    # Wake
    woke_panel  = wake_panel(panels, bound, α)
    woke_vector = wake_vector(woke_panel, panels)
    woke_matrix = [ -woke_vector zeros(length(panels), length(panels) -2) woke_vector ]

    # AIC
    AIC     = doublet_matrix(panels, panels) + woke_matrix
    boco    = dot.(collocation_point.(panels), Ref(u)) - woke_vector * dot(u, r_te)

    # AIC
    # AIC   = influence_matrix(panels, woke_panel)
    # boco  = boundary_vector(ifelse(sources, panels, collocation_point.(panels)), u, r_te) - [ woke_vector; 0 ] .* dot(u, r_te)

    AIC \ boco, AIC, boco
end

"""
    surface_velocities(panels, φs, u, sources :: Bool)

Compute the surface speeds and panel distances given the array of `Panel2D`s, their associated doublet strengths ``φ``s, the velocity ``u``, and a condition whether to disable source terms (``σ = 0``).
"""
function surface_velocities(φs, Δrs, θs, u, sources :: Bool)
    # Δrs   = midpair_map(distance, panels)
    # Δφs   = -midpair_map(-, φs[1:end-1])

    Δφs  = @views φs[1:end-1] - φs[2:end]
    vels = map(zip(Δφs, Δrs, θs)) do (Δφ, Δr, θ)
        if sources
            Δu = dot(u, θ)
        else 
            Δu = 0.
        end
        Δφ / Δr + Δu
    end

    return vels
end

## WAKE VERSIONS
#==========================#

"""
    solve_linear(panels, u, wakes)

Solve the linear aerodynamic system given the array of Panel2Ds, a velocity ``\\vec U``, a vector of wake `Panel2D`s, and an optional named bound for the length of the wake.

The system of equations ``A[φ] = [\\vec{U} ⋅ n̂] - B[σ]`` is solved, where ``A`` is the doublet influence coefficient matrix, ``φ`` is the vector of doublet strengths, ``B`` is the source influence coefficient matrix, and ``σ`` is the vector of source strengths.
"""
function solve_linear(panels :: AbstractArray{<:AbstractPanel2D}, u, wakes)
    AIC  = influence_matrix(panels, wakes)
    boco = boundary_vector(panels, u, [0., 0.])

    AIC \ boco, AIC, boco
end

function solve_linear(panels :: DenseArray{<:AbstractPanel3D}, fs, wakes)
    V∞ = velocity(fs)

    AIC = influence_matrix(panels, wakes)
    boco = boundary_vector(panels, wakes, V∞)

    return AIC \ boco, AIC, boco
end

function influence_matrix(panels :: DenseArray{<:AbstractPanel3D}, wakes)
    panelview = @view permutedims(panels)[:]

    # Foil-Foil interactions
    AIC_ff = doublet_matrix(panelview, panelview)

    # Wake-Foil interactions
    AIC_wf = doublet_matrix(panelview, wakes)

    # Kutta condition
    [   AIC_ff     AIC_wf  ;
      kutta_condition(length(panels), length(wakes)) ]
end

# function influence_matrix!(AIC, sw_panels, wakes)

#     # Foil-foil interaction
#     # for ind in CartesianIndices(sw_panels)
#     #     # i is spanwise index
#     #     i, j = ind.I
#     #     AIC[i,j] = doublet_influence(sw_panels[i], sw_panels[j])
#     # end

#     Ns, Nc = size(sw_panels)

#     # Wake-foil interaction
#     for i in axes(sw_panels, 1)
#         for j in axes(sw_panels, 2)
#             AIC[i,j] = doublet_influence(sw_panels[i], sw_panels[j])
#         end
#         for j in eachindex(wakes)
#             AIC[i,Nc + j] = doublet_influence(sw_panels[i], wakes[j])
#         end
#     end

# end