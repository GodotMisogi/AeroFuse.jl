## Stability derivatives prediction
#==========================================================================================#

# Derivative labels
const Comp = @SLArray (9,7) (# Forces and moment values
        :CX,:CY,:CZ,:Cl,:Cm,:Cn,:CDiff,:CYff,:CLff, 
        # Now their derivatives, using abbreviations for typing convenience
        :CX_M,:CY_M,:CZ_M,:Cl_M,:Cm_M,:Cn_M,:CDiff_M,:CYff_M,:CLff_M, # Mach
        :CX_a,:CY_a,:CZ_a,:Cl_a,:Cm_a,:Cn_a,:CDiff_a,:CYff_a,:CLff_a, # Angle of attack
        :CX_b,:CY_b,:CZ_b,:Cl_b,:Cm_b,:Cn_b,:CDiff_b,:CYff_b,:CLff_b, # Sideslip
        :CX_pb,:CY_pb,:CZ_pb,:Cl_pb,:Cm_pb,:Cn_pb,:CDiff_pb,:CYff_pb,:CLff_pb, # Roll rates 
        :CX_qb,:CY_qb,:CZ_qb,:Cl_qb,:Cm_qb,:Cn_qb,:CDiff_qb,:CYff_qb,:CLff_qb, # Pitch rates
        :CX_rb,:CY_rb,:CZ_rb,:Cl_rb,:Cm_rb,:Cn_rb,:CDiff_rb,:CYff_rb,:CLff_rb  # Yaw rates
    )

# ...and assemble into components
label_derivatives(data, names, Comp = Comp) = @views NamedTuple(names[i] => Comp(data[:,i,:]) for i in Base.axes(data, 2))

# Convert speed and rates into non-dimensional  coefficients
scale_freestream(fs :: Freestream, refs :: References) = 
    [ 
        refs.speed / refs.sound_speed; # Mach number
        fs.alpha;
        fs.beta;
        rate_coefficient(fs, refs)
    ]

# Closure to generate results with input vector
function freestream_derivatives!(y, x, aircraft, fs, ref, name, axes = Wind())
    # Scale Mach number
    ref = @views setproperties(ref, 
        speed = x[1] * ref.sound_speed
    )

    # Transform angles and scale rotation vector
    fs = @views setproperties(fs, 
        alpha = rad2deg(x[2]), 
        beta = rad2deg(x[3]), 
        omega = x[4:end] ./ rate_coefficient(1, ref.speed, ref.span, ref.chord)
    )

    # Solve system
    system = solve_case(
        aircraft, fs, ref,
        name = name
    )

    # Evaluate nearfield coefficients
    CFs, CMs = surface_coefficients(system; axes = axes)

    # Evaluate farfield coefficients
    q = dynamic_pressure(system.reference.density, system.reference.speed)
    FFs = map(farfield_forces(system)) do F
        force_coefficient(F, q, system.reference.area)
    end
    
    # Create array of nearfield and farfield coefficients as a column, for each component as a row.
    comp_coeffs = @views mapreduce(name -> [ sum(CFs[name]); sum(CMs[name]); FFs[name] ], hcat, keys(system.vortices))

    # Massive hack to perform summation for multiple components
    if length(size(comp_coeffs)) > 1
        y[:,1:end-1] = comp_coeffs
        y[:,end] = sum(comp_coeffs, dims = 2)
    else
        y[:,1] = comp_coeffs
        y[:,end] = comp_coeffs
    end

    # return [ comp_coeffs sum_coeffs ] # Append sum of all forces for aircraft
end

function freestream_derivatives(aircraft, fs, ref; axes = Wind(), name = :aircraft, print = false, print_components = false, farfield = false)
    # Reference values and scaling inputs
    x = scale_freestream(fs, ref)

    names     = [ reduce(vcat, keys(aircraft)); name ]
    num_comps = length(names)

    y = zeros(eltype(x), 9, num_comps)
    ∂y∂x = zeros(eltype(x), 9, num_comps, length(x))
    jacobian!(∂y∂x, (y, x) -> freestream_derivatives!(y, x, aircraft, fs, ref, name, axes), y, x)
    
    data = cat(y, ∂y∂x, dims = 3) # Appending outputs with derivatives

    # Labelled array for convenience
    comps = label_derivatives(data, names)

    # Printing
    if print_components
        @views [ print_derivatives(comps[comp], comp, farfield) for comp in keys(comps) ]
    elseif print
        @views print_derivatives(comps[name], name, farfield)
    end

    return comps
end

"""
    freestream_derivatives(
        system :: VortexLatticeSystem,
        name = :aircraft,
        axes = Wind(),
        print = false,
        print_components = false
    )

Obtain the force and moment coefficients of a `VortexLatticeAnalysis` and their derivatives with respect to freestream values.
"""
freestream_derivatives(system :: VortexLatticeSystem; axes = Wind(), name = :aircraft, print = false, print_components = false, farfield = false) = freestream_derivatives(system.vortices, system.freestream, system.reference; axes = axes, name = name, print = print, print_components = print_components, farfield = farfield)

## Linearized stability analysis
#==========================================================================================#

"""
    longitudinal_stability_derivatives(dvs, U, m, Iyy, Q, S, c)

Compute the stability derivatives for the forces and moments in the longitudinal plane ``[X,Z,M]_{u,w,q}``.

The inputs are force and moment coefficient stability derivatives matrix `dvs`, the freestream speed `U`, the mass `m` and longitudinal moment of inertia ``I_{yy}``, the dynamic pressure Q, reference area ``S`` and chord length ``c``.
"""
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

"""
    lateral_stability_derivatives(dvs, U, m, Iyy, Q, S, c)

Compute the stability derivatives for the forces and moments in the lateral plane ``[Y,L,N]_{v,p,r}``.

The inputs are force and moment coefficient stability derivatives matrix `dvs`, the freestream speed `U`, the mass `m` and lateral moment sof inertia ``I_{xx}, I_{zz}``, the dynamic pressure Q, reference area ``S`` and span length ``b``.
"""
function lateral_stability_derivatives(dvs, U, m, Ixx, Izz, Q, S, b)
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