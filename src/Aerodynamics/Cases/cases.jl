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
    solve_case(horseshoe_panels :: Matrix{Panel3D}, normals, freestream :: Freestream, r_ref, ρ = 1.225; symmetry = false)

Evaluate a vortex lattice case given an array of `Panel3D`s with associated normal vectors (not necessarily the same as the panels' normals), a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, with an option for symmetry.
"""
function solve_case(horseshoe_panels :: Array{<: Panel3D}, normals, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref, r_ref)
    # Make horseshoes and collocation points
    horseshoes = horseshoe_line.(horseshoe_panels)

    # Solve system
    Γs = reshape(solve_system(horseshoes[:], normals[:], U, Ω), size(horseshoe_panels))

    # Compute forces and moments
    surface_forces, surface_moments, trefftz_force = case_dynamics(Γs, horseshoes, U, α, β, Ω, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(surface_forces, surface_moments, trefftz_force, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref)

    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs
end

"""
    solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref; span_num :: Integer = 5, chord_num :: Integer = 10)

Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
"""
function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, r_ref = [0.25, 0, 0], area_ref = projected_area(wing), chord_ref = mean_aerodynamic_chord(wing), span_ref = span(wing), mu_ref = 1.5e-5, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, viscous = false, a_ref = 330., x_tr = 0.3, spacing = spanwise_spacing(wing))
    # Unpack Freestream
    U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Determine spanwise panel distribution and spacing
    space     = ifelse(typeof(spacing) <: String, [spacing], spacing)
    span_nums = number_of_spanwise_panels(wing, span_num)

    # Compute panels and normals
    horseshoe_panels, camber_panels = vlmesh_wing(wing, span_nums, chord_num, space)
    normals = panel_normal.(camber_panels)

    # Compute forces and moments
    nearfield_coeffs, farfield_coeffs, CFs, CMs, horseshoes, Γs = solve_case(horseshoe_panels, normals, U, α, β, Ω, rho_ref, area_ref, chord_ref, span_ref, r_ref)
    
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

function solve_case(components :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft", print = false, print_components = false) where T <: Real
    # Get panels
    meshes = values(components)
    
    # Flattening for VLM
    horseshoe_panels = first.(meshes)
    normals          = last.(meshes)
    horsies          = reduce(vcat, vec.(horseshoe_panels))
    normies          = reduce(vcat, vec.(normals))

    # Get required vortex lattice variables, i.e. horseshoes, collocation points and normals
    horseshoes = horseshoe_line.(horsies)
    horseshoes_arr = [ horseshoe_line.(horses) for horses in horseshoe_panels ]
    
    # Unpack Freestream
    U, α, β, Ω = aircraft_velocity(freestream), freestream.alpha, freestream.beta, freestream.omega

    # Solve system
    Γs = solve_system(horseshoes, normies, U, Ω)

    # Reshaping
    panel_sizes = size.(horseshoe_panels)
    panel_inds 	= [ 0; cumsum(prod.(panel_sizes)) ]
    Γs_arr 		= reshape_array(Γs, panel_inds, panel_sizes)

    # Compute forces and moments
    results = case_dynamics.(Γs_arr, horseshoes_arr, Ref(Γs), Ref(horseshoes), Ref(U), α, β, Ref(Ω), rho_ref, Ref(r_ref))
    forces  = getindex.(results, 1)
    moments = getindex.(results, 2) 
    trefftz = getindex.(results, 3)

    # Components' non-dimensional forces and moments
    data = evaluate_coefficients.(forces, moments, trefftz, Ref(U), α, β, Ref(Ω), rho_ref, area_ref, chord_ref, span_ref)

    nf_comp_coeffs = getindex.(data, 1)
    ff_comp_coeffs = getindex.(data, 2)
    CFs            = getindex.(data, 3)
    CMs            = getindex.(data, 4)

    # Aircraft's non-dimensional forces and moments
    nf_coeffs = reduce((x, y) -> x .+ y, nf_comp_coeffs) # Sum nearfield coefficients
    ff_coeffs = reduce((x, y) -> x .+ y, ff_comp_coeffs) # Sum farfield coefficients
    name_CFs  = reduce(vcat, vec.(CFs))                  # Collect surface force coefficients
    name_CMs  = reduce(vcat, vec.(CMs))                  # Collect surface moment coefficients

    # Printing
    if print_components; print_coefficients.(nf_comp_coeffs, ff_comp_coeffs, keys(components)) end
    if print; 			 print_coefficients(nf_coeffs, ff_coeffs, name) 					   end

    # Dictionary assembly
    name_data = (nf_coeffs, ff_coeffs, name_CFs, name_CMs, horsies, normies, horseshoes, Γs)
    comp_data = tuple.(nf_comp_coeffs, ff_comp_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes_arr, Γs_arr)

    names 	= [ name						 ; # Aircraft name
                (collect ∘ keys)(components) ] # Component names
    data  	= [ name_data ;	# Aircraft data
                comp_data ]	# Component data

    Dict(names .=> data)
end

# Mutating version
function solve_case(aircraft :: Dict{String, Tuple{Matrix{Panel3D{T}}, Matrix{SVector{3,T}}}}, state :: VLMState) where T <: Real
    # Build surfaces and systems
    system, surfs = build_system(aircraft)

    # Solve case
    solve_case!(system, surfs, state)

    system, surfs
end

## Method extensions from submodules
#==========================================================================================#

streamlines(freestream :: Freestream, points, horseshoes, Γs, length, num_steps :: Integer) = VortexLattice.streamlines.(points, Ref(velocity(freestream)), Ref(freestream.omega), Ref(horseshoes), Ref(Γs), Ref(length), Ref(num_steps))

VortexLattice.VLMState(fs :: Freestream{<: Real}; r_ref = zeros(3), rho_ref = 1.225, area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") = VLMState(fs.V, fs.alpha, fs.beta, fs.omega, r_ref = r_ref, rho_ref = rho_ref, area_ref = area_ref, chord_ref = chord_ref, span_ref = span_ref, name = name)

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

function print_coefficients(surfs, state :: VLMState{T}) where T <: Real
    coeffs = aerodynamic_coefficients(surfs, state)
    print_coefficients.(first.(values(coeffs)), last.(values(coeffs)), (collect ∘ keys)(coeffs))
end

function print_coefficients(surf :: VLMSurface{T}, state :: VLMState{T}) where T <: Real
    nf, ff = aerodynamic_coefficients(surf, state)
    print_coefficients(nf, ff, name(surf))
end