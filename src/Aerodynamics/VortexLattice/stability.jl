## Stability derivatives prediction
#==========================================================================================#

# Derivative labels
const Derivs = @SLArray (9,7) (# Forces and moment values
        :CX,:CY,:CZ,:Cl,:Cm,:Cn,:CDi,:CY,:CL, 
        # Now their derivatives, using abbreviations for typing convenience
        :CX_Ma,:CY_Ma,:CZ_Ma,:Cl_Ma,:Cm_Ma,:Cn_Ma,:CDi_Ma,:CY_Ma,:CL_Ma, # Mach number
        :CX_al,:CY_al,:CZ_al,:Cl_al,:Cm_al,:Cn_al,:CDi_al,:CY_al,:CL_al, # Angle of attack
        :CX_be,:CY_be,:CZ_be,:Cl_be,:Cm_be,:Cn_be,:CDi_be,:CY_be,:CL_be, # Sideslip
        :CX_pb,:CY_pb,:CZ_pb,:Cl_pb,:Cm_pb,:Cn_pb,:CDi_pb,:CY_pb,:CL_pb, # Roll rates 
        :CX_qb,:CY_qb,:CZ_qb,:Cl_qb,:Cm_qb,:Cn_qb,:CDi_qb,:CY_qb,:CL_qb, # Pitch rates
        :CX_rb,:CY_rb,:CZ_rb,:Cl_rb,:Cm_rb,:Cn_rb,:CDi_rb,:CY_rb,:CL_rb  # Yaw rates
    )

# ...and assemble into components
label_derivatives(data, names, Derivs = Derivs) = @views NamedTuple(names[i] => Derivs(data[:,i,:]) for i in Base.axes(data, 2))

# Convert speed and rates into non-dimensional coefficients
scale_freestream(fs :: Freestream, refs :: References) = [ 
        mach_number(refs); # Mach number
        fs.alpha;
        fs.beta;
        rate_coefficient(# Non-dimensionalize rates
            geometry_to_stability_axes(# Transform from geometry to stability axes
                flip_xz(fs.omega), # Flip signs of x- and z- components
                fs.alpha # Angle of rotation (attack)
            ), 
            refs.speed, refs.span, refs.chord # Non-dimensional factors
        )
    ]

# Closure to generate results with input vector
@views function freestream_derivatives!(y, x, aircraft, fs, ref, compressible, axes)
    # Scale Mach number
    ref = @views setproperties(ref, 
        speed = x[1] * ref.sound_speed
    )

    # Define new freestream with scaled and transformed
    fs = @views setproperties(fs, 
        alpha = x[2], 
        beta = x[3], 
        omega = flip_xz(# Flip signs of x- and z- components
            stability_to_geometry_axes(# Transform from stability to geometry axes
                x[4:end] ./ rate_coefficient(1, ref.speed, ref.span, ref.chord), # Dimensionalize rates
                x[2] # Angle of rotation (attack)
            )
        )
    )

    # Solve system
    system = VortexLatticeSystem(aircraft, fs, ref, compressible)

    # Evaluate nearfield coefficients
    CFs, CMs = surface_coefficients(system; axes)

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

    return nothing
end

function freestream_derivatives(aircraft, fs, ref; axes = Stability(), name = :aircraft, compressible = false, print = false, print_components = false, farfield = false)
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
        @views [ print_derivatives(comps[comp], comp; farfield) for comp in keys(comps) ]
    elseif print
        @views print_derivatives(comps[name], name; farfield)
    end

    return comps
end

"""
    freestream_derivatives(
        system :: VortexLatticeSystem;
        axes :: AbstractAxisSystem = Stability(),
        name = :aircraft,
        print = false,
        print_components = false,
        farfield = false
    )

Compute the force and moment coefficients of the components of a `VortexLatticeSystem` and their derivatives with respect to freestream values: Mach ``M`` (if compressible), angles of attack ``α`` and sideslip ``β``, and non-dimensionalized rotation rates ``p̄, q̄, r̄`` (in stability axes).

The axes of the force and moment coefficients can be changed by passing any `AbstractAxisSystem` (such as `Body(), Geometry(), Wind(), Stability()`) to the named `axes` argument. The nearfield force and moment coefficients are reported in stability axes by default. Note that the derivatives with respect to the rotation rates will still refer to the rotation vector in stability axes ``(p̄, q̄, r̄)``.

Optional printing arguments are provided for the components and the entire system, along with the corresponding farfield coefficients if needed.
"""
freestream_derivatives(system :: VortexLatticeSystem; axes = Stability(), name = :aircraft, print = false, print_components = false, farfield = false) = freestream_derivatives(system.vortices, system.freestream, system.reference; compressible = system.compressible, axes, name, print, print_components, farfield)