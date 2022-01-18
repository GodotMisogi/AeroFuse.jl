## Printing
#==========================================================================================#

"""
    print_info(wing :: AbstractWing)

Print the relevant geometric characteristics of a `HalfWing` or `Wing`.
"""
function print_info(wing :: AbstractWing, head = "")
    labels = [ "Aspect Ratio", "Span (m)", "Area (m²)", "Mean Aerodynamic Chord (m)", "Mean Aerodynamic Center (m)" ]
    wing_info = properties(wing)
    data = Any[ labels wing_info ]
    header = [ head, "Value" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true)

    pretty_table(data, header, alignment = [:c, :c], highlighters = h1, vlines = :none, formatters = ft_round(8))
end

"""
    print_coefficients(nf_coeffs, ff_coeffs, name = "")

Print a pretty table of the nearfield and farfield coefficients with an optional name.
"""
function print_coefficients(nf_coeffs :: AbstractVector{T}, ff_coeffs :: AbstractVector{T}, name = "") where T <: Real
    coeffs = [ ifelse(length(nf_coeffs) == 8, ["CD", "CDv"], []); [ "CDi", "CY", "CL", "Cl", "Cm", "Cn" ] ]
    data = [ coeffs nf_coeffs [ ff_coeffs; fill("—", 3) ] ]
    head = [ name, "Nearfield", "Farfield" ]
    h1 = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true)

    pretty_table(data, head, alignment = [:c, :c, :c], highlighters = h1, vlines = :none, formatters = ft_round(8))
end

"""
    print_derivatives(nf_coeffs, ff_coeffs, name = ""; 
                      axes = "")

Print a pretty table of the aerodynamic coefficients and derivatives with an optional name and axis name.
"""
function print_derivatives(comp, name = ""; axes = "")
    coeffs  = ["CX", "CY", "CZ", "Cℓ", "Cm", "Cn"]
    nf_vars = [ "$name" "Values" "" "" "Aerodynamic" "Derivatives" "" "" ; "" "$axes" "∂/∂U, m⁻¹s" "∂/∂α, 1/rad" "∂/∂β, 1/rad" "∂/∂p̄" "∂/∂q̄" "∂/∂r̄" ]
    nf_rows = @views [ coeffs comp[1:6,:] ]

    pretty_table(nf_rows, nf_vars, alignment = :c, header_crayon = Crayon(bold = true), subheader_crayon = Crayon(foreground = :yellow, bold = true), highlighters = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true), vlines = :none, formatters = ft_round(8))
end

"""
    print_coefficients(system :: VLMSystem, name = :aircraft;
                       components = false)

Print a pretty table of the total nearfield and farfield coefficients of a `VLMSystem` with an optional name.

A named Boolean argument `components` is provided to alsoenable the printing of any possible components.
"""
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