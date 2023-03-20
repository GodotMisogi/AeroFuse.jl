## Stability derivatives prediction
#==========================================================================================#

# Derivative labels
const Derivs = @SLArray (9,7) (# Forces and moment values
        :CX,:CY,:CZ,:Cl,:Cm,:Cn,:CDiff,:CYff,:CLff, 
        # Now their derivatives, using abbreviations for typing convenience
        :CX_Ma,:CY_Ma,:CZ_Ma,:Cl_Ma,:Cm_Ma,:Cn_Ma,:CDiff_Ma,:CYff_Ma,:CLff_Ma, # Mach number
        :CX_al,:CY_al,:CZ_al,:Cl_al,:Cm_al,:Cn_al,:CDiff_al,:CYff_al,:CLff_al, # Angle of attack
        :CX_be,:CY_be,:CZ_be,:Cl_be,:Cm_be,:Cn_be,:CDiff_be,:CYff_be,:CLff_be, # Sideslip
        :CX_pb,:CY_pb,:CZ_pb,:Cl_pb,:Cm_pb,:Cn_pb,:CDiff_pb,:CYff_pb,:CLff_pb, # Roll rates 
        :CX_qb,:CY_qb,:CZ_qb,:Cl_qb,:Cm_qb,:Cn_qb,:CDiff_qb,:CYff_qb,:CLff_qb, # Pitch rates
        :CX_rb,:CY_rb,:CZ_rb,:Cl_rb,:Cm_rb,:Cn_rb,:CDiff_rb,:CYff_rb,:CLff_rb  # Yaw rates
    )

# ...and assemble into components
label_derivatives(data, names, Derivs = Derivs) = @views NamedTuple(names[i] => Derivs(data[:,i,:]) for i in Base.axes(data, 2))

# Convert speed and rates into non-dimensional coefficients
scale_freestream(fs :: Freestream, refs :: References) = 
    [ 
        mach_number(refs); # Mach number
        fs.alpha;
        fs.beta;
        rate_coefficient(fs, refs)
    ]

# Closure to generate results with input vector
function freestream_derivatives!(y, x, aircraft, fs, ref, compressible = false, axes = Wind())
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
    system = VortexLatticeSystem(aircraft, fs, ref, compressible)

    # Evaluate nearfield coefficients
    CFs, CMs = surface_coefficients(system; axes = axes)

    # Evaluate farfield coefficients
    FFs = farfield_coefficients(system)

    # Create array of nearfield and farfield coefficients as a column, for each component as a row.
    comp_coeffs = @views mapreduce(name -> [ sum(CFs[name]); sum(CMs[name]); FFs[name] ], hcat, keys(system.vortices))

    # Massive hack to perform summation of coefficients from multiple components
    if length(size(comp_coeffs)) > 1
        y[:,1:end-1] = comp_coeffs
        y[:,end] = sum(comp_coeffs, dims = 2)
    else
        y[:,1] = comp_coeffs
        y[:,end] = comp_coeffs
    end

    # return nothing
end

function freestream_derivatives(aircraft, fs, ref; axes = Wind(), name = :aircraft, compressible = false, print = false, print_components = false, farfield = false)
    # Reference values and scaling inputs
    x = scale_freestream(fs, ref)

    names = [ reduce(vcat, keys(aircraft)); name ]
    num_comps = length(names)

    y = zeros(eltype(x), 9, num_comps)
    ∂y∂x = zeros(eltype(x), 9, num_comps, length(x))
    jacobian!(∂y∂x, (y, x) -> freestream_derivatives!(y, x, aircraft, fs, ref, compressible, axes), y, x)
    
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

Obtain the force and moment coefficients of a `VortexLatticeSystem` and the derivatives of its components with respect to freestream values.
"""
freestream_derivatives(system :: VortexLatticeSystem; axes = Wind(), name = :aircraft, print = false, print_components = false, farfield = false) = freestream_derivatives(system.vortices, system.freestream, system.reference; axes = axes, name = name, compressible = system.compressible, print = print, print_components = print_components, farfield = farfield)