function solve_case(horseshoe_panels, normals, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225; symmetry = false)
    # Solve system
    @timeit "Solve System" Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], normals[:], freestream, symmetry), size(horseshoe_panels)...)

    geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment = case_dynamics(Γs, horseshoes, freestream, r_ref, ρ, symmetry)

    geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, Γs, horseshoes
end

function case_coefficients(wing :: Union{Wing, HalfWing}, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, V, ρ)
    nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, V, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, V, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    nearfield_coeffs, farfield_coeffs
end

symmetric_case_coefficients(wing :: Union{Wing, HalfWing}, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, V, ρ) = case_coefficients(wing, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, V, 2ρ)

function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, ρ :: Real, r_ref = [0.25, 0., 0.]; span_num :: Integer = 5, chord_num :: Integer = 10)
    if typeof(wing) == Wing && wing.left === wing.right && freestream.β == 0. && freestream.Ω == zeros(3)
        # Compute panels and normals
        @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing.right, span_num, chord_num)

        normals = panel_normal.(camber_panels)

        geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, Γs, horseshoes = solve_case(horseshoe_panels, normals, freestream, r_ref, ρ, symmetry = true)

        nearfield_coeffs, farfield_coeffs = symmetric_case_coefficients(wing.right, force, moment, drag / 2, trans_rates, trefftz_force, trefftz_moment, freestream.mag, ρ)

        Γs = [ Γs[:,end:-1:1] Γs ]
        horseshoe_panels = [ reflect_xz.(horseshoe_panels[:,end:-1:1]) horseshoe_panels ]
        camber_panels = [ reflect_xz.(camber_panels[:,end:-1:1]) camber_panels ]
        horseshoes = [ reflect_xz.(horseshoes[:,end:-1:1]) horseshoes ]
    else
        # Compute panels and normals
        @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
        normals = panel_normal.(camber_panels)
        
        geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, Γs, horseshoes = solve_case(horseshoe_panels, normals, freestream, r_ref, ρ)

        nearfield_coeffs, farfield_coeffs = case_coefficients(wing, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, freestream.mag, ρ)
    end
    
    nearfield_coeffs, farfield_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs
end

function solve_case(foil :: Foil, freestream :: Uniform2D; num_panels :: Integer = 60)
    # @timeit "Make Panels" 
    airfoil = paneller(foil, num_panels)
    
    # @timeit "Solve Case" 
    φs, cl = solve_problem(airfoil, Laplace.velocity(freestream))
    
    cl
end


# function solve_case(components :: Dict{<: Aircraft, Tuple{Integer, Integer}}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225)
#     # Compute panels
#     @timeit "Meshing" meshes = [ vlmesh_wing(comp, size_panels...) for (comp, size_panels) in components ] 
#     @timeit "Make Panels" horseshoe_panels, camber_panels = first.(meshes), last.(meshes)
    
#     # Solve system with normalised velocities
#     @timeit "Solve System" Γs, horseshoes = solve_horseshoes(horseshoe_panels[:], camber_panels[:], freestream, symmetry)

#     Γs = reshape(Γs, size(horseshoe_panels)...)

#     # Compute nearfield dynamics
#     @timeit "Nearfield Dynamics" geom_forces, geom_moments = nearfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)
#     force, moment = sum(geom_forces), sum(geom_moments)
#     drag = nearfield_drag(force, freestream)

#     @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels[:]))

#     @timeit "Transforming Axes" trans_forces, trans_moments, trans_rates = body_to_wind_axes(force, moment, freestream.Ω, freestream)

#     # Compute farfield dynamics
#     @timeit "Farfield Dynamics" trefftz_force, trefftz_moment = farfield_dynamics(Γs, horseshoes, freestream, r_ref, ρ)

#     # Compute non-dimensional coefficients
#     @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

#     @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, freestream.mag, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

#     nearfield_coeffs, farfield_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs
# end
