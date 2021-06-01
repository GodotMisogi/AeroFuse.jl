## Stability derivative cases
function scale_inputs(freestream :: Freestream, b, c)
	# Building input vector
	x0 = [ rad2deg(freestream.α); 
		   rad2deg(freestream.β); 
		   rate_coefficient(freestream.Ω, freestream.V, b, c) ]

	# Unscaling non-dimensional rate coefficients
	scale = 2freestream.V .* [1/b, 1/c, 1/b]

	x0, scale
end

function print_case(data, comp)
	nf, ff, derivs = data[comp]
	print_coefficients(comp, nf, ff)
	print_derivatives(comp, derivs)
end

function solve_stability_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, span_num :: Union{Integer, Vector{<: Integer}}, chord_num :: Integer, name = "Wing", viscous = false, x_tr = 0.3, print = false)
	# Reference values and scaling inputs
	S, b, c = area_ref, span_ref, chord_ref
	x, scale = scale_inputs(freestream, b , c)

	# Closure to generate results with input vector
	function stab(x)
		fs 	 = Freestream(freestream.V, x[1], x[2], x[3:end] .* scale)
        data = solve_case(wing, fs;
						  rho_ref   = rho_ref, 
						  r_ref     = r_ref, 
						  area_ref  = S, 
						  span_ref  = b, 
						  chord_ref = c, 
						  chord_num = chord_num,
						  span_num  = span_num,
						  viscous 	= viscous,
						  x_tr 		= x_tr)

        [ data[1]; data[2] ] 
    end

	y = ifelse(viscous, zeros(16), zeros(12))
	result = DiffResults.JacobianResult(y, x)
	result = ForwardDiff.jacobian!(result, stab, x)

	vars 	= DiffResults.value(result)
	derivs 	= DiffResults.jacobian(result)
	nf 		= ifelse(viscous, vars[1:11], vars[1:9])
	ff 		= ifelse(viscous, vars[12:end], vars[10:end])
	dvs 	= ifelse(viscous, derivs[3:end,:], derivs)

	if print
		print_coefficients(name, nf, ff)
		print_derivatives(name, dvs)
	end

	nf, ff, dvs
end

function solve_stability_case(aircraft :: AbstractDict{String, NTuple{2, Matrix{Panel3D{T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft", print = false, print_components = false) where T <: Real
	# Reference values and scaling inputs
	S, b, c = area_ref, span_ref, chord_ref
	x, scale = scale_inputs(freestream, b , c)

	# Closure to generate results with input vector
    function stab(x)
		fs 	 = Freestream(freestream.V, x[1], x[2], x[3:end] .* scale)
        data = solve_case(aircraft, fs,
            			  rho_ref   = rho_ref, 
            			  r_ref     = r_ref, 
            			  area_ref  = S, 
            			  span_ref  = b, 
            			  chord_ref = c,
						  name 		= name)

		# Creates array of nearfield and farfield coefficients for each component as a row vector.
		coeffs = hcat((let comps = data[name]; [ comps[1]; comps[2] ] end for name in names)...)
    end

	names 	= [ keys(aircraft)..., name ]
	num_comps = length(names)

	y 		= zeros(12, num_comps)
	result 	= DiffResults.JacobianResult(y, x)
	result 	= ForwardDiff.jacobian!(result, stab, x)

	vars 	= DiffResults.value(result)
	derivs 	= DiffResults.jacobian(result)

	# Dictionary assembly
	ranges  = (1:num_comps) .* 12
	bounds  = zip([ 1; ranges[1:end-1] .+ 1], ranges)
	data 	= Dict(name => (vars[1:9, i], vars[10:end, i], derivs[first(inds):last(inds),:]) for (i, (name, inds)) in (enumerate ∘ zip)(names, bounds)) 
	
	if print
		print_case(data, name)
	end
	
	if print_components
		print_case.(Ref(data), names)
	end

	data
end


