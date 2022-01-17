## Stability derivative cases
function scale_inputs(fs :: Freestream, refs :: References)
    # Scaling rate coefficients
    scale = [refs.span, refs.chord, refs.span] ./ (2 * refs.speed)

    # Building input vector
    x0 = [ refs.speed;
           fs.alpha;
           fs.beta;
           fs.omega .* scale ]

    x0, scale
end

"""
    solve_case_derivatives(components :: Vector{Horseshoe}, fs :: Freestream, refs :: References;
                           name = :aircraft :: Symbol, 
                           print = false :: Boolean,
                           print_components = false :: Boolean)

Obtain the values and derivatives of a vortex lattice analysis given a vector of `Horseshoe`s, a `Freestream` condition, and `Reference` values.
"""
function solve_case_derivatives(aircraft, fs :: Freestream, ref :: References; name = :aircraft, print = false, print_components = false)
    # Reference values and scaling inputs
    x, scale = scale_inputs(fs, ref)

    # Closure to generate results with input vector
    function freestream_derivatives(x)
        # @set ref.speed = x[1]
        ref_c  = setproperties(ref, speed = x[1])
        fs_c   = setproperties(fs, alpha = rad2deg(x[2]), beta = rad2deg(x[3]), omega = x[4:end] ./ scale)
        system = solve_case(aircraft, fs_c, ref_c,
                            name      = name)

        NFs = nearfield_coefficients(system)
        FFs = farfield_coefficients(system)

        # Creates array of nearfield and farfield coefficients for each component as a row vector.
        comp_coeffs = mapreduce(name -> [ NFs[name]; FFs[name] ], hcat, valkeys(system.vortices))
        all_coeffs  = [ comp_coeffs sum(comp_coeffs, dims = 2) ] # Append sum of all forces for aircraft
    end

    names     = [ reduce(vcat, keys(aircraft)); name ]
    num_comps = length(names)

    y       = zeros(9, num_comps)
    result  = JacobianResult(y, x)
    jacobian!(result, freestream_derivatives, x)

    vars    = value(result)
    derivs  = jacobian(result)

    # Reshaping
    data  = cat(vars, reshape(derivs, 9, num_comps, 6), dims = 3)

    # Labelled array for convenience
    Comp = @SLArray (9,7) (:CD,:CY,:CL,:Cl,:Cm,:Cn,:CD_ff,:CY_ff,:CL_ff,:CD_speed,:CY_speed,:CL_speed,:Cl_speed,:Cm_speed,:Cn_speed,:CD_ff_speed,:CY_ff_speed,:CL_ff_speed,:CD_alpha,:CY_alpha,:CL_alpha,:Cl_alpha,:Cm_alpha,:Cn_alpha,:CD_ff_alpha,:CY_ff_alpha,:CL_ff_alpha,:CD_beta,:CY_beta,:CL_beta,:Cl_beta,:Cm_beta,:Cn_beta,:CD_ff_beta,:CY_ff_beta,:CL_ff_beta,:CD_pb,:CY_pbar,:CL_pbar,:Cl_pbar,:Cm_pbar,:Cn_pbar,:CD_ff_pbar,:CY_ff_pbar,:CL_ff_pbar,:CD_qbar,:CY_qbar,:CL_qbar,:Cl_qbar,:Cm_qbar,:Cn_qbar,:CD_ff_qbar,:CY_ff_qbar,:CL_ff_qbar,:CD_rbar,:CY_rbar,:CL_rbar,:Cl_rbar,:Cm_rbar,:Cn_rbar,:CD_ff_rbar,:CY_ff_rbar,:CL_ff_rbar)

    comps = @views NamedTuple(names[i] => Comp(data[:,i,:]) for i in axes(data, 2))

    # Printing
    if print_components
        @views [ print_derivatives(comps[comp], comp) for comp in keys(comps) ]
    elseif print
        @views print_derivatives(comps[name], name)
    end

    comps
end

## Printing
#==========================================================================================#

function print_case(nf, ff, derivs, comp)
    print_coefficients(nf, ff, comp)
    print_derivatives(derivs, comp)
end

function print_case(data, comp)
    nf, ff, derivs = data[comp]
    print_case(nf, ff, derivs, comp)
end

function print_derivatives(comp, name = ""; browser = false)
    coeffs  = ["CD", "CY", "CL", "Cl", "Cm", "Cn", "CD_ff", "CY_ff", "CL_ff"]
    nf_vars = [ "$name" "Values" "" "" "Aerodynamic" "Derivatives" "" "" ; "" "" "∂/∂U, m⁻¹s" "∂/∂α, 1/rad" "∂/∂β, 1/rad" "∂/∂p̄" "∂/∂q̄" "∂/∂r̄" ]
    nf_rows = [ coeffs comp ]

    if browser
        pretty_table(HTML, nf_rows, nf_vars, alignment = :c, tf = tf_html_minimalist, backend = :html, highlighters = HTMLHighlighter( (data,i,j) -> (j == 1), HTMLDecoration(color = "blue", font_weight = "bold")), formatters = ft_round(8))
    else
        pretty_table(nf_rows, nf_vars, alignment = :c, tf = tf_compact, header_crayon = Crayon(bold = true), subheader_crayon = Crayon(foreground = :yellow, bold = true), highlighters = Highlighter( (data,i,j) -> (j == 1), foreground = :cyan, bold = true), vlines = :none, formatters = ft_round(8))
    end
end