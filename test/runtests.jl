using AeroMDAO
using Test

@testset "Kulfan CST Doublet-Source Panel Method" begin
    alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
    alpha_l = [-0.2, -0.1, -0.1, -0.001]
    dzs     = (0., 0.)

    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
    
    uniform = Uniform2D(1., 5.)
    cl, cls, cms, cps, panels = solve_case(airfoil, uniform; num_panels = 80)

    @test isapprox(cl, 0.84188988, atol = 1e-6)
    @test isapprox(sum(cls), 0.84073703, atol = 1e-6)
    @test isapprox(sum(cms), -0.26104277, atol = 1e-6)
end

@testset "Airfoil Processing" begin
    # Import and read airfoil coordinates
    foilpath = joinpath((dirname ∘ dirname ∘ pathof)(AeroMDAO), "test/CRM.dat")
    coords   = read_foil(foilpath)

    # Cosine spacing
    cos_foil = cosine_foil(coords, 51)

    # Split airfoil
    up, low  = split_foil(cos_foil)

    # Convert coordinates to Kulfan CST variables
    num_dv   = 4
    alpha_u, alpha_l = coords_to_CST(up, num_dv), coords_to_CST(low, num_dv)

    # Generate same airfoil using Kulfan CST parametrisation
    cst_foil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (0., 0.), 0.0)

    uniform  = Uniform2D(1., 5.)
    cl, cls, cms, cps, panels = solve_case(cst_foil, uniform; num_panels = 80)

    @test isapprox(cl, 0.85736965, atol = 1e-6)
    @test isapprox(sum(cls), 0.85976886, atol = 1e-6)
    @test isapprox(sum(cms), -0.29766116, atol = 1e-6)
end

@testset "NACA-4 Wing Vortex Lattice Method" begin
    foil        =   naca4((2,4,1,2))
    wing_right  =   HalfWing(Foil.(foil for i ∈ 1:3),
                            [0.18, 0.16, 0.08],
                            [2., 0., -2.],
                            [0.5, 0.5],
                            [0., 11.3],
                            [1.14, 8.])

    wing = Wing(wing_right, wing_right)

    ρ = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω = [0.0, 0.0, 0.0]
    uniform = Freestream(10.0, 2.0, 2.0, Ω)
    nf_coeffs, ff_coeffs, dv_coeffs = solve_stability_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5)

    display(dv_coeffs)

    CL_nf, CDi_nf, CY_nf, Cl_nf, Cm_nf, Cn_nf, p_b_nf, q_b_nf, r_b_nf = nf_coeffs
    CL_ff, CDi_ff, CY_ff = ff_coeffs

    nf_tests = [0.00094437, -0.00049066, 0.11137544, -0.001868206, -0.00440614, -9.807e-5, 0.0, 0.0, 0.0]
    ff_tests = [0.00103027, -0.00053694, 0.11136510]
    dv_tests = [ 0.000405506   5.54895e-7    0.00233895  -0.0108045     0.000566369 ;
                 5.35821e-5   -0.000242975   0.157886     0.000211956  -0.019372600 ;
                 0.02732       1.21838e-5   -0.0340142    0.414198      0.005331970 ;
                 0.000157786  -0.000924136   0.737322     0.00279386   -0.087240700 ;
                -0.000564098  -6.71773e-5    0.0267506   -0.0132855    -0.003225710 ;
                -3.14918e-5   -4.92097e-5    0.0306778   -0.000785897   0.002498710 ] 

    [ @test isapprox(nf_c, nf_t, atol = 1e-6) for (nf_c, nf_t) in zip(nf_coeffs, nf_tests) ]
    [ @test isapprox(ff_c, ff_t, atol = 1e-6) for (ff_c, ff_t) in zip(ff_coeffs, ff_tests) ]
    [ @test isapprox(dv_c, dv_t, atol = 1e-6) for (dv_c, dv_t) in zip(dv_coeffs[1:6,:], dv_tests) ]
end