## Pretty-printing and show methods
#==========================================================================================#

# Axis systems
Base.show(io :: IO, :: Geometry)  = print(io, "Geometry")
Base.show(io :: IO, :: Body)      = print(io, "Body")
Base.show(io :: IO, :: Stability) = print(io, "Stability")
Base.show(io :: IO, :: Wind)      = print(io, "Wind")

# Reference values
function Base.show(io :: IO, refs :: References)
    println(io, "References: ")
    for fname in fieldnames(typeof(refs))
        println(io, "    ", fname, " = ", getfield(refs, fname))
    end
end

function Base.show(io :: IO, ring :: AbstractVortex)
    println(io, "Vortex: ")
    for fname in fieldnames(typeof(ring))
        println(io, "    ", fname, " = ", getfield(ring, fname))
    end
end

# Vortex lattice system
function Base.show(io :: IO, sys :: VortexLatticeSystem)     
    println(io, "VortexLatticeSystem -")
    println(io, length(sys.vortices), " ", eltype(sys.vortices), " Elements")
    println(io, "Axes: ", sys.axes)
    show(io, sys.freestream)
    show(io, sys.reference)
end


"""
    print_coefficients(nf_coeffs, ff_coeffs, name = "")

Print a pretty table of the nearfield and farfield coefficients with an optional name.
"""
function print_coefficients(nf_coeffs, ff_coeffs, name = "")
    visc = ifelse(length(nf_coeffs) == 8, ["CD", "CDv"], [])
    coeffs = [ [ visc; "CX"; "CY"; "CZ"; "Cl"; "Cm"; "Cn" ] [ visc; "CDi"; "CYff"; "CL"; ""; ""; "" ] ]
    data = [ coeffs[:,1] [ nf_coeffs... ] coeffs[:,2] [ [ ff_coeffs...]; fill("—", 3) ] ]
    head = [ name, "Nearfield", "", "Farfield" ]
    h1 = Highlighter( (data,i,j) -> (j == 1) || (j == 3), foreground = :cyan, bold = true)

    pretty_table(data, 
        header = head, 
        alignment = [:c, :c, :c, :c], 
        highlighters = h1, 
        vlines = :none, 
        formatters = ft_round(8)
    )
end

"""
    print_derivatives(nf_coeffs, ff_coeffs, name = ""; 
                      axes = "")

Print a pretty table of the aerodynamic coefficients and derivatives with an optional name and axis name.
"""
function print_derivatives(comp, name = "", farfield = false; axes = "")
    coeffs  = ["CX", "CY", "CZ", "Cℓ", "Cm", "Cn", "CDi ff", "CY ff", "CL ff"]
    nf_vars = (["$name" "Values" "" "" "Freestream" "Derivatives" "" ""], ["" "$axes" "∂/∂M" "∂/∂α, 1/rad" "∂/∂β, 1/rad" "∂/∂p̄" "∂/∂q̄" "∂/∂r̄" ])
    ff_index = ifelse(farfield, 9, 6)
    nf_rows = @views [ coeffs[1:ff_index] comp[1:ff_index,:] ]

    pretty_table(nf_rows, 
        header = nf_vars, 
        alignment = :c, 
        header_crayon = Crayon(bold = true), 
        subheader_crayon = Crayon(foreground = :yellow, bold = true), 
        highlighters = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true), 
        vlines = :none, formatters = ft_round(8)
    )
end

"""
    print_coefficients(system :: VortexLatticeSystem, name = :aircraft;
                       components = false)

Print a pretty table of the total nearfield and farfield coefficients of a `VortexLatticeSystem` with an optional name.

A named Boolean argument `components` is provided to also enable the printing of any possible components.
"""
function print_coefficients(system :: VortexLatticeSystem, name = :aircraft; components = false)
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