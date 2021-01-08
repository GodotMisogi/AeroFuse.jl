# function solve_case(horseshoe_panels :: AbstractVector{Panel3D}, camber_panels :: AbstractVector{Panel3D}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225; symmetry = false)

# end

function case_coefficients(wing :: Union{Wing, HalfWing}, force :: SVector{3,<: Real}, moment :: SVector{3,<: Real}, drag :: Real, trans_rates :: SVector{3,<: Real}, trefftz_force :: SVector{3,<: Real}, trefftz_moment :: SVector{3,<: Real}, V :: Real, ρ :: Real)
    @timeit "Nearfield Coefficients" nearfield_coeffs = aerodynamic_coefficients(force, moment, drag, trans_rates, V, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    @timeit "Farfield Coefficients" farfield_coeffs = aerodynamic_coefficients(trefftz_force, trefftz_moment, trans_rates, V, projected_area(wing), span(wing), mean_aerodynamic_chord(wing), ρ)

    nearfield_coeffs, farfield_coeffs
end

symmetric_case_coefficients(wing :: Union{Wing, HalfWing}, force :: SVector{3,<: Real}, moment :: SVector{3,<: Real}, drag :: Real, trans_rates :: SVector{3,<: Real}, trefftz_force :: SVector{3,<: Real}, trefftz_moment :: SVector{3,<: Real}, V :: Real, ρ :: Real) = case_coefficients(wing, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, V, 2ρ)

function solve_symmetric_case(wing :: HalfWing, freestream :: Freestream, span_num :: Integer, chord_num :: Integer, r_ref, ρ)
    # Compute panels
    @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
    
    # Solve system
    @timeit "Solve System" Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], camber_panels[:], freestream, true), size(horseshoe_panels)...)

    geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment = case_dynamics(Γs, horseshoes, freestream, r_ref, ρ, true)
    
    # @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels[:]))

    nearfield_coeffs, farfield_coeffs = symmetric_case_coefficients(wing, force, moment, drag / 2, trans_rates, trefftz_force, trefftz_moment, freestream.mag, ρ)

    nearfield_coeffs, farfield_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs
end

function solve_case(wing :: Union{Wing, HalfWing}, freestream :: Freestream, r_ref = SVector(0.25, 0., 0.), ρ = 1.225; span_num :: Integer = 15, chord_num :: Integer = 5)
    if wing.left === wing.right && freestream.β == 0.
        solve_symmetric_case(wing.right, freestream, span_num, chord_num, r_ref, ρ)
    else
        # Compute panels
        @timeit "Make Panels" horseshoe_panels, camber_panels = vlmesh_wing(wing, span_num, chord_num)
        
        # Solve system
        @timeit "Solve System" Γs, horseshoes = reshape.(solve_horseshoes(horseshoe_panels[:], camber_panels[:], freestream, false), size(horseshoe_panels)...)

        geom_forces, geom_moments, force, moment, drag, trans_rates, trefftz_force, trefftz_moment = case_dynamics(Γs, horseshoes, freestream, r_ref, ρ)
        
        # @timeit "Pressure Distribution" cps = pressure_coefficient.(geom_forces, ρ, freestream.mag, panel_area.(camber_panels[:]))

        nearfield_coeffs, farfield_coeffs = case_coefficients(wing, force, moment, drag, trans_rates, trefftz_force, trefftz_moment, freestream.mag, ρ)

        nearfield_coeffs, farfield_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs
    end
end

function solve_case(foil :: Foil, freestream :: Uniform2D, num_panels :: Integer = 60)
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
