function solve_case(horseshoe_panels :: Array{Panel3D}, camber_panels :: Array{Panel3D}, freestream :: Freestream, Ω = SVector(0., 0., 0.), r_ref = SVector(0.25, 0., 0.), ρ = 1.225; print = true, symmetry = false)
    vel = velocity(freestream)

    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry)

    # Compute near-field forces
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, VortexLattice.area.(camber_panels))
    # println(cps)

    force, moment = sum(geom_forces), sum(geom_moments)
    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, Ω, freestream)
    drag = nearfield_drag(force, freestream)

    # @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    # # Compute non-dimensional coefficients
    # @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    force, drag, moment, horseshoes, Γs
    # horseshoe_panels, camber_panels, horseshoes, Γs
    # coeffs
end

function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, Ω = SVector(0., 0., 0.), r_ref = SVector(0.25, 0., 0.), ρ = 1.225; span_num :: Integer = 15, chord_num :: Integer = 5, print = true)
    vel = velocity(freestream)

    # Experimental: Symmetry condition
    symmetry = false
    # if typeof(wing) == Wing
    #     println("Symmetric Case")
    #     symmetry, wing = wing.left === wing.right ? (true, wing.right) : (false, wing)
    # end

    # reset_timer!()

    # Compute panels
    @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
    
    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry)

    # Compute near-field forces
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels))
    # println(cps)

    force, moment = sum(geom_forces), sum(geom_moments)
    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream.Ω, freestream)
    drag = nearfield_drag(force, freestream)

    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    # Compute non-dimensional coefficients
    @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    # # Optionally print data
    @timeit "Printing" if print 
        println("Nearfield:") 
        print_dynamics(nearfield_coeffs...)
        println("\nFarfield:")
        print_dynamics(farfield_coeffs...)
    end

    horseshoe_panels, camber_panels, horseshoes, Γs
    # nearfield_coeffs
end

# solve_case(wings :: Array{Aircraft}, freestream :: Freestream, r_ref = (0.25, 0, 0), ρ = 1.225; span_num = 15, chord_num = 5, print = true) =  

function solve_case(panels :: Array{Panel2D}, freestream :: Uniform2D)
    # Compute doublet and source strengths
    φs, σs = solve_strengths(panels, freestream)
    freestream_speed = norm(velocity(freestream))
    pressure_coeffs = pressure_coefficient.(freestream_speed, panel_velocities(panels, freestream, φs[1:end-1]))

    # Make panels for plotting
    # dub_src_panels = DoubletSourcePanel2D.(panels, φs[1:end-1], σs, pressure_coeffs)

    # Compute lift coefficient
    diff_pans = [ panel_dist(panel_1, panel_2) for (panel_1, panel_2) ∈ (collect ∘ eachrow ∘ midgrad)(panels) ]
    cl = sum(lift_coefficient.(pressure_coeffs, diff_pans ./ 2, panel_angle.(panels)))

    cl
end
