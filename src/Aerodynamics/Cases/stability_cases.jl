## Stability derivatives prediction
#==========================================================================================#

function scale_inputs(fs :: Freestream, refs :: References)
    # Scaling rate coefficients
    scale = (refs.span, refs.chord, refs.span) ./ (2 * refs.speed)

    # Building input vector
    x0 = [ refs.speed;
           fs.alpha;
           fs.beta;
           fs.omega .* scale ]

    x0, scale
end

"""
    solve_case_derivatives(components :: Vector{Horseshoe}, 
                           fs :: Freestream, 
                           refs :: References;
                           name = :aircraft, 
                           print = false,
                           print_components = false)

Obtain the values and derivatives of a vortex lattice analysis given a vector of `Horseshoe`s, a `Freestream` condition, and `Reference` values.
"""
function solve_case_derivatives(aircraft, fs :: Freestream, ref :: References; axes = Wind(), name = :aircraft, print = false, print_components = false)
    # Reference values and scaling inputs
    x, scale = scale_inputs(fs, ref)

    # Closure to generate results with input vector
    function freestream_derivatives(x)
        ref_c  = setproperties(ref, speed = x[1])
        fs_c   = setproperties(fs, alpha = rad2deg(x[2]), beta = rad2deg(x[3]), omega = x[4:end] ./ scale)
        system = solve_case(aircraft, fs_c, ref_c,
                            name      = name)

        CFs, CMs = surface_coefficients(system; axes = axes)
        FFs = farfield_coefficients(system)
        
        # Create array of nearfield and farfield coefficients for each component as a row vector.
        comp_coeffs = mapreduce(name -> [ sum(CFs[name]); sum(CMs[name]); FFs[name] ], hcat, keys(system.vortices))
        
        [ comp_coeffs sum(comp_coeffs, dims = 2) ] # Append sum of all forces for aircraft
    end

    names     = [ reduce(vcat, keys(aircraft)); name ]
    num_comps = length(names)

    y       = zeros(9, num_comps)
    result  = JacobianResult(y, x)
    jacobian!(result, freestream_derivatives, x)

    values = value(result)
    derivs = jacobian(result)
    data   = cat(values, reshape(derivs, 9, num_comps, 6), dims = 3) # Reshaping

    # Labelled array for convenience
    Comp = @SLArray (9,7) (:CX,:CY,:CZ,:Cl,:Cm,:Cn,:CD_ff,:CY_ff,:CL_ff,:CX_speed,:CY_speed,:CZ_speed,:Cl_speed,:Cm_speed,:Cn_speed,:CD_ff_speed,:CY_ff_speed,:CL_ff_speed,:CX_alpha,:CY_alpha,:CZ_alpha,:Cl_alpha,:Cm_alpha,:Cn_alpha,:CD_ff_alpha,:CY_ff_alpha,:CL_ff_alpha,:CX_beta,:CY_beta,:CZ_beta,:Cl_beta,:Cm_beta,:Cn_beta,:CD_ff_beta,:CY_ff_beta,:CL_ff_beta,:CX_pb,:CY_pbar,:CZ_pbar,:Cl_pbar,:Cm_pbar,:Cn_pbar,:CD_ff_pbar,:CY_ff_pbar,:CL_ff_pbar,:CX_qbar,:CY_qbar,:CZ_qbar,:Cl_qbar,:Cm_qbar,:Cn_qbar,:CD_ff_qbar,:CY_ff_qbar,:CL_ff_qbar,:CX_rbar,:CY_rbar,:CZ_rbar,:Cl_rbar,:Cm_rbar,:Cn_rbar,:CD_ff_rbar,:CY_ff_rbar,:CL_ff_rbar)

    comps = @views NamedTuple(names[i] => Comp(data[:,i,:]) for i in Base.axes(data, 2))

    # Printing
    if print_components
        @views [ print_derivatives(comps[comp], comp) for comp in keys(comps) ]
    elseif print
        @views print_derivatives(comps[name], name)
    end

    comps
end

## Linearized stability analysis
#==========================================================================================#

function longitudinal_stability_derivatives(dvs, U, m, Iyy, Q, S, c)
    QS  = Q * S
    c1  = QS / (m * U) 
    c2  = QS / (Iyy * U)

    X_u = c1 * dvs.CX_speed
    Z_u = c1 * dvs.CZ_speed
    M_u = c2 * dvs.Cm_speed * c

    X_w = c1 * dvs.CX_alpha
    Z_w = c1 * dvs.CZ_alpha
    M_w = c2 * dvs.Cm_alpha * c

    X_q = c1 / 2 * dvs.CX_qbar
    Z_q = c1 / 2 * dvs.CZ_qbar
    M_q = c2 / 2 * dvs.Cm_qbar * c

    X_u, Z_u, M_u, X_w, Z_w, M_w, X_q, Z_q, M_q
end

longitudinal_stability_matrix(X_u, Z_u, M_u, X_w, Z_w, M_w, X_q, Z_q, M_q, U₀, g) = 
    [ X_u X_w        0         -g
      Z_u Z_w     U₀ + Z_q      0
      M_u M_w       M_q         0
       0   0         1          0 ]

function lateral_stability_derivatives(dvs, U, m, Iyy, Izz, Q, S, b)
    QS  = Q * S
    c1  = QS / (m * U) 
    c2  = QS / (Ixx * U)
    c3  = QS / (Izz * U)

    Y_v = c1 * dvs.CY_beta
    L_v = c2 * dvs.Cl_beta * b
    N_v = c3 * dvs.Cn_beta * b

    Y_p = c1 / 2 * dvs.CY_pbar * b
    L_p = c2 / 2 * dvs.Cl_pbar * b^2
    N_p = c3 / 2 * dvs.Cn_pbar * b^2

    Y_r = c1 / 2 * dvs.CY_rbar * b
    L_r = c2 / 2 * dvs.Cl_rbar * b^2
    N_r = c3 / 2 * dvs.Cn_rbar * b^2

    Y_v, L_v, N_v, Y_p, L_p, N_p, Y_r, L_r, N_r
end

lateral_stability_matrix(Y_v, L_v, N_v, Y_p, L_p, N_p, Y_r, L_r, N_r, U₀, θ₀, g) =
    [ Y_v Y_p (Y_r - U₀) (g * cos(θ₀))
      L_v L_p     L_r         0      
      N_v N_p     N_r         0      
       0   1       0          0       ]