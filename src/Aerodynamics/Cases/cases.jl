## 3D cases
#==========================================================================================#

"""
    print_info(wing :: AbstractWing)

Print the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
function print_info(wing :: AbstractWing, head = ""; browser = false)
    labels = [ "Aspect Ratio", "Span (m)", "Area (m²)", "Mean Aerodynamic Chord (m)", "Mean Aerodynamic Center (m)" ]
    wing_info = properties(wing)
    data = Any[ labels wing_info ]
    header = [ head, "Value" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true)
    if browser
        pretty_table(String, data, header, alignment = [:c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, header, alignment = [:c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end

"""
    solve_case(components :: Vector{Horseshoe}, fs :: Freestream, refs :: References;
               name = :aircraft :: Symbol, 
               print = false :: Boolean,
               print_components = false :: Boolean)

Perform a vortex lattice analysis given a vector of `Horseshoe`s, a `Freestream` condition, and `Reference` values.
"""
function solve_case(components, freestream :: Freestream, refs :: References; name = :aircraft, print = false, print_components = false, finite_core = false)
    system = solve_system(components, freestream, refs, finite_core)

    # Printing if needed
    if print_components
        print_coefficients(system, name, components = true)
    elseif print
        print_coefficients(system, name)
    end

    system
end

function print_coefficients(system :: VLMSystem, name = :aircraft; components = false)
    if components
        nf_c = nearfield_coefficients(system)
        ff_c = farfield_coefficients(system) 
        [ print_coefficients(nf_c[key], ff_c[key], key) for key in keys(system.vortices) ]
        print_coefficients(nearfield(system), farfield(system), name)
    else
        print_coefficients(nearfield(system), farfield(system), name)
    end

    nothing
end

## Placeholder for functions I'm not sure where to put
#==========================================================================================#

## Span-loading
function span_loads(panels, CFs_wind, S)
    CFs  = combinedimsview(CFs_wind)
    CDis = @views CFs[1,:,:]
    CYs  = @views CFs[2,:,:]
    CLs  = @views CFs[3,:,:]

    area_scale = S ./ vec(sum(panel_area, panels, dims = 1))
    span_CDis  = vec(sum(CDis, dims = 1)) .* area_scale
    span_CYs   = vec(sum(CYs,  dims = 1)) .* area_scale
    span_CLs   = vec(sum(CLs,  dims = 1)) .* area_scale

    ys = @views vec(mean(x -> midpoint(x)[2], panels, dims = 1))

    [ ys span_CDis span_CYs span_CLs ]
end

span_loads(wing :: WingMesh, CFs_wind, S) = span_loads(chord_panels(wing), CFs_wind, S)

## Mesh connectivities
triangle_connectivities(inds) = @views [ 
                                         vec(inds[1:end-1,1:end-1]) vec(inds[1:end-1,2:end]) vec(inds[2:end,2:end])    ;
                                         vec(inds[2:end,2:end])     vec(inds[2:end,1:end-1]) vec(inds[1:end-1,1:end-1])
                                       ]

## Extrapolating surface values to neighbouring points
function extrapolate_point_mesh(mesh, weight = 0.75)
    m, n   = size(mesh)
    points = zeros(eltype(mesh), m + 1, n + 1)

    # The quantities are measured at the bound leg (0.25×)
    @views points[1:end-1,1:end-1] += weight       / 2 * mesh
    @views points[1:end-1,2:end]   += weight       / 2 * mesh
    @views points[2:end,1:end-1]   += (1 - weight) / 2 * mesh
    @views points[2:end,2:end]     += (1 - weight) / 2 * mesh

    points
end

## Printing
#==========================================================================================#

function print_coefficients(nf_coeffs :: AbstractVector{T}, ff_coeffs :: AbstractVector{T}, name = ""; browser = false) where T <: Real
    coeffs = [ ifelse(length(nf_coeffs) == 8, ["CD", "CDp"], []); [ "CDi", "CY", "CL", "Cl", "Cm", "Cn" ] ]
    data = [ coeffs nf_coeffs [ ff_coeffs; fill("—", 3) ] ]
    head = [ name, "Nearfield", "Farfield" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true)
    if browser
        pretty_table(HTML, data, head, alignment = [:c, :c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, head, alignment = [:c, :c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end


# """
#     solve_case(wing :: AbstractWing, freestream :: Freestream, ρ :: Real, r_ref; span_num :: Integer = 5, chord_num :: Integer = 10)

# Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``ρ`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
# """
# function solve_case(wing :: AbstractWing, freestream :: Freestream, refs = References(); mu_ref = 1.5e-5, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, viscous = false, a_ref = 330., x_tr = 0.3, spacing = symmetric_spacing(wing))
#     # Determine spanwise panel distribution and spacing
#     span_nums = number_of_spanwise_panels(wing, span_num)

#     # Compute panels and normals
#     horseshoe_panels = mesh_horseshoes(wing, span_nums, chord_num; spacings = spacing)
#     camber_panels    = mesh_cambers(wing, span_nums, chord_num; spacings = spacing)
#     normals          = panel_normal.(camber_panels)

#     # Make horseshoes and collocation points
#     components = ComponentVector(wing = Horseshoe.(horseshoe_panels, normals))

#     # Evaluate case
#     system = solve_system(components, body_frame_velocity(freestream), freestream.omega, finite_core)

#     # Viscous drag evaluation
#     if viscous
#         # Compute profile drag via wetted-area equivalent skin-friction method
#         CDv = profile_drag_coefficient(wing, x_tr, refs.speed, refs.sound_speed, refs.area, refs.dynamic_viscosity)

#         # Compute profile drag using local dissipation based on power-balance method of Sato.
#         # u_es = surface_velocities(system)
#         # CDv = local_dissipation_drag(wing, panel_area.(camber_panels), ρ_es, u_es, x_tr, V, ρ, M, μ)

#         # Add profile and total drag coefficients
#         # nf_coeffs = [ nearfield_coeffs[1] + CDv; CDv; nearfield_coeffs ]
#         # ff_coeffs = [  farfield_coeffs[1] + CDv; CDv; farfield_coeffs  ]
#     else
#         # nf_coeffs = nearfield_coeffs
#         # ff_coeffs = farfield_coeffs
#     end

#     system
# end