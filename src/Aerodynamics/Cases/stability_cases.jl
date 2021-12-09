## Stability derivative cases
function scale_inputs(freestream :: Freestream, b, c)
    # Building input vector
    x0 = [ freestream.alpha;
           freestream.beta;
           rate_coefficient(freestream.omega, freestream.V, b, c) ]

    # Unscaling non-dimensional rate coefficients
    scale = 2freestream.V .* [1/b, 1/c, 1/b]

    x0, scale
end

function print_case(nf, ff, derivs, comp)
    print_coefficients(nf, ff, comp)
    print_derivatives(derivs, comp)
end

function print_case(data, comp)
    nf, ff, derivs = data[comp]
    print_case(nf, ff, derivs, comp)
end

function solve_stability_case(wing :: AbstractWing, freestream :: Freestream; rho_ref = 1.225, area_ref = projected_area(wing), chord_ref = mean_aerodynamic_chord(wing), r_ref = [ 0.25 * chord_ref, 0., 0.], span_ref = span(wing), span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, name = "Wing", viscous = false, x_tr = 0.3, print = false, spacing = symmetric_spacing(wing))
    # Reference values and scaling inputs
    S, b, c = area_ref, span_ref, chord_ref
    x, scale = scale_inputs(freestream, b, c)

    # Closure to generate results with input vector
    function stab(x)
        fs   = Freestream(freestream.V, rad2deg(x[1]), rad2deg(x[2]), x[3:end] .* scale)
        data = solve_case(wing, fs;
                          rho_ref   = rho_ref,
                          r_ref     = r_ref,
                          area_ref  = S,
                          span_ref  = b,
                          chord_ref = c,
                          chord_num = chord_num,
                          span_num  = span_num,
                          viscous   = viscous,
                          x_tr      = x_tr,
                          spacing   = spacing)

        [ data[1]; data[2] ] # Concatenate nearfield and farfield coefficients for DiffResults value
    end

    # Set up Jacobian system and evaluate
    y = ifelse(viscous, zeros(13), zeros(9))
    result = DiffResults.JacobianResult(y, x)
    result = ForwardDiff.jacobian!(result, stab, x)
    vars   = DiffResults.value(result)
    derivs = DiffResults.jacobian(result)

    # Gather results
    nf  = ifelse(viscous, vars[1:8], vars[1:6])       # Nearfield coefficients
    ff  = ifelse(viscous, vars[9:end], vars[7:end])  # Farfield coefficients
    dvs = ifelse(viscous, derivs[[1,4:8...],:], derivs[1:6,:])  # Nearfield stability derivatives, uses total drag for viscous cases just in case of generalisations

    # Print if needed
    if print; print_case(nf, ff, derivs, name) end

    nf, ff, dvs
end

function solve_stability_case(aircraft, freestream :: Freestream, ref :: References; name = :aircraft, print = false, print_components = false)
    # Reference values and scaling inputs
    b, c = ref.span, ref.chord
    x, scale = scale_inputs(freestream, b , c)

    # Closure to generate results with input vector
    function stab(x)
        fs   = Freestream(freestream.V, rad2deg(x[1]), rad2deg(x[2]), x[3:end] .* scale)
        data = solve_case(aircraft, fs, ref,
                          name      = name)

        NFs = nearfield_coefficients(data)
        FFs = farfield_coefficients(data)

        # Creates array of nearfield and farfield coefficients for each component as a row vector.
        comp_coeffs = mapreduce(name -> [ NFs[name]; FFs[name] ], hcat, valkeys(data.horseshoes))
        all_coeffs  = [ comp_coeffs sum(comp_coeffs, dims = 2) ] # Append sum of all forces for aircraft
    end

    names     = [ reduce(vcat, keys(aircraft)); name ]
    num_comps = length(names)

    y       = zeros(9, num_comps)
    result  = DiffResults.JacobianResult(y, x)
    result  = ForwardDiff.jacobian!(result, stab, x)

    vars    = DiffResults.value(result)
    derivs  = DiffResults.jacobian(result)

    # Reshaping
    data  = cat(vars, reshape(derivs, 9, num_comps, 5), dims = 3)
    comps = @views NamedTuple(names[i] => ComponentArray(NF  = data[1:6,i,1],
                                                         dNF = data[1:6,i,2:end], 
                                                         FF  = data[7:end,i,1], 
                                                         dFF = data[7:end,i,2:end]) 
                              for i in axes(data, 2))

    # Printing
    if print_components
        @views [ print_derivatives(comps[comp], comp) for comp in keys(comps) ]
    elseif print
        @views print_derivatives(comps[name], name)
    end

    comps
end


