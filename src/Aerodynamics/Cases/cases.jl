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
function solve_case(horseshoe_panels :: Array{<: Panel3D}, normals, freestream :: Freestream, ρ, r_ref = zeros(3))
    U = aircraft_velocity(freestream)
    Ω = freestream.Ω

    # Solve system
    Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], normals[:], U, Ω), size(horseshoe_panels)...)

    # Compute forces and moments
    geom_forces, geom_moments, trefftz_force = case_dynamics(Γs, horseshoes, U, freestream.α, freestream.β, Ω, ρ, r_ref)
    
    geom_forces, geom_moments, trefftz_force, Γs, horseshoes
end

"""
    solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref; span_num :: Integer = 5, chord_num :: Integer = 10)

Evaluate a vortex lattice case given a `Wing` or `HalfWing` with a `Freestream`, reference density ``\\rho`` and reference point ``r_\\text{ref}`` for moments, ``n_s`` span-wise panels and ``n_c`` chord-wise panels.
"""
function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, r_ref = [0.25, 0, 0], area_ref = 1, chord_ref = 1, span_ref = 1, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, viscous = false, a_ref = 330., x_tr = 0.3)
    # Determine spanwise panel distribution
    if typeof(span_num) <: Integer
        span_nums = ifelse(typeof(wing) <: Wing, ceil.(Int, span_num / 2 .* wing.right.spans / span(wing.right)), span_num)
    else
        span_nums = ceil.(Int, span_num / 2)
    end

    # Compute panels and normals
    horseshoe_panels, camber_panels = vlmesh_wing(wing, span_nums, chord_num)
    normals = panel_normal.(camber_panels)

    # Compute forces and moments
    panel_forces, panel_moments, trefftz_force, Γs, horseshoes = solve_case(horseshoe_panels, normals, freestream, rho_ref, r_ref)

    # Compute aerodynamic coefficients
    nearfield_coeffs, farfield_coeffs, CFs, CMs = evaluate_coefficients(panel_forces, panel_moments, trefftz_force, aircraft_velocity(freestream), freestream.α, freestream.β, freestream.Ω, rho_ref, area_ref, chord_ref, span_ref)
    
    if viscous
        # Compute profile drag via wetted-area equivalent skin-friction method
        wing_right = ifelse(typeof(wing) <: Wing, wing.right, wing)
        μ = 1.5e-5
        Dp_by_q = wetted_area_drag(wing_right, x_tr, freestream.V, rho_ref, a_ref, μ)

        CDv = ifelse(typeof(wing) <: Wing, 2Dp_by_q, Dp_by_q) / area_ref

        # Add profile drag for printing
        nf_coeffs = [ nearfield_coeffs[1] + CDv; CDv; nearfield_coeffs ]
        ff_coeffs = [ farfield_coeffs[1] + CDv; CDv; farfield_coeffs ]
    else
        nf_coeffs = nearfield_coeffs
        ff_coeffs = farfield_coeffs
    end

    nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes, Γs
end

function solve_case(components :: Dict{String, NTuple{2, Matrix{Panel3D{T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft", print = false, print_components = false) where T <: Real
    # Get panels
    meshes = values(components)
    
    # Flattening for VLM
    horseshoe_panels = first.(meshes)
    camber_panels 	 = last.(meshes)

    cammies = vcat(vec.(camber_panels)...)
    horsies = vcat(vec.(horseshoe_panels)...)
    normals = panel_normal.(cammies)
    
    # Solve system
    U = aircraft_velocity(freestream)
    Ω = freestream.Ω
       Γs, horseshoes = solve_horseshoes(horsies, normals, U, Ω)

    # Reshaping
    panel_sizes = size.(horseshoe_panels)
    panel_inds 	= [ 0; cumsum(prod.(panel_sizes)) ]
    
    Γs_arr 		   = reshape_array(Γs, panel_inds, panel_sizes)
    horseshoes_arr = reshape_array(horseshoes, panel_inds, panel_sizes)

    # Compute forces and moments
    results = case_dynamics.(Γs_arr, horseshoes_arr, Ref(Γs), Ref(horseshoes), Ref(U), freestream.α, freestream.β, Ref(Ω), rho_ref, Ref(r_ref))

    # Components' non-dimensional forces and moments
    data = [ evaluate_coefficients(dyn..., U, freestream.α, freestream.β, Ω, rho_ref, area_ref, chord_ref, span_ref) for dyn in results ]

    nf_comp_coeffs 	= getindex.(data, 1)
    ff_comp_coeffs 	= getindex.(data, 2)
    CFs 			= getindex.(data, 3)
    CMs 			= getindex.(data, 4)

    # Aircraft's non-dimensional forces and moments
    nf_coeffs = reduce((x, y) -> x .+ y, nf_comp_coeffs)
    ff_coeffs = reduce((x, y) -> x .+ y, ff_comp_coeffs)
    name_CFs  = vcat(vec.(CFs)...)
    name_CMs  = vcat(vec.(CMs)...)

    # Printing
    if print_components
        print_coefficients.(keys(components), nf_comp_coeffs, ff_comp_coeffs)
    end

    if print
        print_coefficients(name, nf_coeffs, ff_coeffs)
    end

    # Dictionary assembly
    name_data = (nf_coeffs, ff_coeffs, name_CFs, name_CMs, horsies, cammies, horseshoes, Γs)
    comp_data = tuple.(nf_comp_coeffs, ff_comp_coeffs, CFs, CMs, horseshoe_panels, camber_panels, horseshoes_arr, Γs_arr)

    names 	= [ name						 ; 
                (collect ∘ keys)(components) ]
    data  	= [ name_data ;
                comp_data ]

    Dict(names .=> data)
end

streamlines(freestream :: Freestream, points, horseshoes, Γs, length, num_steps :: Integer) = VortexLattice.streamlines.(points, Ref(velocity(freestream)), Ref(freestream.Ω), Ref(horseshoes), Ref(Γs), Ref(length), Ref(num_steps))