## 3D cases
#==========================================================================================#

"""
    print_info(wing :: Union{Wing, HalfWing})

Print the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
function print_info(wing :: Union{Wing, HalfWing}, head = ""; browser = false)
    labels = [ "Span (m)", "Area (m²)", "MAC (m)", "Aspect Ratio"]
    wing_info = [ info(wing)... ]
    data = [ labels wing_info ]
    header = [ head, "Value" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true)
    if browser
        pretty_table(String, data, header, alignment = [:c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, header, alignment = [:c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end

"""
    solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref; span_num :: Integer = 5, chord_num :: Integer = 10)

Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
"""
function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, area_ref = projected_area(wing), chord_ref = mean_aerodynamic_chord(wing), r_ref = [0.25 * chord_ref, 0., 0.], span_ref = span(wing), mu_ref = 1.5e-5, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, viscous = false, a_ref = 330., x_tr = 0.3, spacing = symmetric_spacing(wing))
    # Unpack Freestream
    U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Determine spanwise panel distribution and spacing
    space     = ifelse(typeof(spacing) <: String, [spacing], spacing)
    span_nums = number_of_spanwise_panels(wing, span_num)

    # Compute panels and normals
    horseshoe_panels = mesh_horseshoes(wing, span_nums, chord_num; spacings = space)
    camber_panels    = mesh_cambers(wing, span_nums, chord_num; spacings = space)
    normals          = panel_normal.(camber_panels)

    # Make horseshoes and collocation points
    horseshoes = Horseshoe.(horseshoe_panels, normals)

    # Compute forces and moments
    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs = evaluate_case(horseshoes, U, α, β, Ω, rho_ref, r_ref, area_ref, chord_ref, span_ref)

    # Viscous drag evaluation
    if viscous
        # Compute profile drag via wetted-area equivalent skin-friction method
        CDv = profile_drag_coefficient(wing, x_tr, freestream.V, rho_ref, a_ref, area_ref, mu_ref)

        # Add profile and total drag coefficients
        nf_coeffs = [ nearfield_coeffs[1] + CDv; CDv; nearfield_coeffs ]
        ff_coeffs = [  farfield_coeffs[1] + CDv; CDv; farfield_coeffs  ]
    else
        nf_coeffs = nearfield_coeffs
        ff_coeffs = farfield_coeffs
    end

    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs
end

# rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1,

function solve_case(components, freestream :: Freestream, refs :: References; name = :aircraft, print = false, print_components = false)
    # Unpack Freestream
    # U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Evaluate case
    results = solve_system(components, aircraft_velocity(freestream), freestream.omega)

    system = VLMSystem(components, results, freestream, refs)

    
    # NamedTuple{keys(components)}(tuple_comp)

    # Printing if needed
    if print_components
        nf_c = nearfield_coefficients(system)
        ff_c = farfield_coefficients(system) 
        print_coefficients(nearfield(system), farfield(system), name)
        [ print_coefficients(nf_c[key], ff_c[key], key) for key in keys(components) ]
    elseif print
        print_coefficients(nearfield(system), farfield(system), name)
    end

    system
end

## State cases

# solve_case(horseshoe_panels :: Array{<: Panel3D}, normals, state :: VLMState) = solve_case(horseshoe_panels, normals, state.speed, state.alpha, state.beta, state.omega, rho_ref = state.rho_ref, r_ref = state.r_ref, area_ref = state.area_ref, chord_ref = state.chord_ref, span_ref = state.span_ref)

# solve_case(components :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, state :: VLMState) = evaluate_case(components, state.speed, state.alpha, state.beta, state.omega, rho_ref = state.rho_ref, r_ref = state.r_ref, area_ref = state.area_ref, chord_ref = state.chord_ref, span_ref = state.span_ref, name)

# Mutating version
# function solve_case(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, state :: VLMState) where T <: Real
#     # Build surfaces and systems
#     system = build_system(aircraft)

#     # Solve case for given state
#     evaluate_case!(system, state)

#     system
# end

## Method extensions from submodules
#==========================================================================================#

streamlines(freestream :: Freestream, points, horseshoes, Γs, length, num_steps :: Integer) = VortexLattice.streamlines.(points, Ref(velocity(freestream)), Ref(freestream.omega), Ref(horseshoes), Ref(Γs), Ref(length), Ref(num_steps))

# VortexLattice.VLMState(fs :: Freestream{<: Real}; r_ref = zeros(3), rho_ref = 1.225, area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") = VLMState(fs.V, fs.alpha, fs.beta, fs.omega, r_ref = r_ref, rho_ref = rho_ref, area_ref = area_ref, chord_ref = chord_ref, span_ref = span_ref, name = name)

## Printing
#==========================================================================================#

function print_coefficients(nf_coeffs :: AbstractVector{T}, ff_coeffs :: AbstractVector{T}, name = ""; browser = false) where T <: Real
    coeffs = [ ifelse(length(nf_coeffs) == 8, ["CD", "CDp"], []); [ "CDi", "CY", "CL", "Cl", "Cm", "Cn" ] ]
    data = [ coeffs nf_coeffs [ ff_coeffs; fill("—", 3) ] ]
    head = [ name, "Nearfield", "Farfield" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true)
    if browser
        pretty_table(String, data, head, alignment = [:c, :c, :c], tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(data, head, alignment = [:c, :c, :c], tf = tf_compact, highlighters = h1, vlines = :none, formatters = ft_round(8))
    end
end

function print_derivatives(derivs, name = ""; browser = false)
    coeffs = ["∂CD", "∂CY", "∂CL", "∂Cl", "∂Cm", "∂Cn"]
    nf_vars = [ "$name" "" "Nearfield" "Stability" "Derivatives" "" ; "" "∂α, 1/rad" "∂β, 1/rad" "∂p̄" "∂q̄" "∂r̄" ]
    nf_rows = [ coeffs derivs ]

    if browser
        pretty_table(String, nf_rows, nf_vars, alignment = :c, tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(color = "blue", font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(nf_rows, nf_vars, alignment = :c, tf = tf_compact, header_crayon = Crayon(bold = true), subheader_crayon = Crayon(foreground = :yellow, bold = true), highlighters = Highlighter( (data,i,j) -> (j == 1), foreground = :blue, bold = true), vlines = :none, formatters = ft_round(8))
    end
end

# function print_coefficients(surfs, state :: VLMState{T}) where T <: Real
#     coeffs = aerodynamic_coefficients(surfs, state)
#     print_coefficients.(first.(values(coeffs)), last.(values(coeffs)), (collect ∘ keys)(coeffs))
# end

# function print_coefficients(surf :: VLMSurface{T}, state :: VLMState{T}) where T <: Real
#     nf, ff = aerodynamic_coefficients(surf, state)
#     print_coefficients(nf, ff, name(surf))
# end