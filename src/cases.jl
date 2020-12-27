function solve_case(horseshoe_panels :: AbstractVector{Panel3D}, camber_panels :: AbstractVector{Panel3D}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225; symmetry = false)
    vel = VortexLattice.velocity(freestream)

    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels, camber_panels, freestream, symmetry)

    # Compute near-field forces
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    force, moment = sum(geom_forces), sum(geom_moments)
    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream)
    drag = nearfield_drag(force, freestream)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels))

    # Compute farfield forces
    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    force, drag, moment, horseshoes, Γs
end

function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225; span_num :: Integer = 15, chord_num :: Integer = 5)
    vel = VortexLattice.velocity(freestream)

    # Experimental: Symmetry condition
    symmetry = false
    # if typeof(wing) == Wing
    #     println("Symmetric Case - \n")
    #     symmetry, wing = wing.left === wing.right ? (true, wing.right) : (false, wing)
    # end

    # Compute panels
    @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
    
    # Solve system
    @timeit "Solve System" Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], camber_panels[:], freestream, symmetry), size(horseshoe_panels)...)

    # Compute near-field dynamics
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs[:], horseshoes[:], freestream, r_ref, ρ)
    force, moment = sum(geom_forces), sum(geom_moments)
    drag = nearfield_drag(force, freestream)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels[:]))

    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream)

    # Compute farfield dynamics
    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    # Compute non-dimensional coefficients
    @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    nearfield_coeffs, farfield_coeffs, cps, horseshoe_panels, camber_panels, horseshoes, Γs
end

function solve_case(components :: Dict{<: Aircraft, Tuple{Integer, Integer}}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225)
    vel = VortexLattice.velocity(freestream)
    # Compute panels
    @timeit "Meshing" meshes = [ vlmesh_wing(comp, size_panels...) for (comp, size_panels) in components ] 
    @timeit "Make Panels" horseshoe_panels, camber_panels = first.(meshes), last.(meshes)
    
    # Solve system with normalised velocities
    @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels[:], camber_panels[:], freestream, symmetry)

    # Compute nearfield dynamics
    @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)
    force, moment = sum(geom_forces), sum(geom_moments)
    drag = nearfield_drag(force, freestream)

    @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels[:]))

    @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream.Ω, freestream)

    # Compute farfield dynamics
    @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

    # Compute non-dimensional coefficients
    @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    nearfield_coeffs, farfield_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs
end

function solve_case(foil :: Foil, freestream :: Uniform2D, num_panels :: Integer = 60)
    @timeit "Make Panels" airfoil = paneller(foil, num_panels)
    
    @timeit "Solve System" φs = solve_strengths(airfoil, freestream)
    
    @timeit "Lift Coefficient" cl = lift_coefficient(airfoil, freestream, φs)
end