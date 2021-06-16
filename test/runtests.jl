using AeroMDAO
using Test

@testset "Kulfan CST Doublet-Source Panel Method" begin
    # Define Kulfan CST coefficients
    alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
    alpha_l = [-0.2, -0.1, -0.1, -0.001]
    dzs     = (0., 0.)

    # Define airfoil
    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
    
    # Define uniform flow
    uniform = Uniform2D(1., 5.)

    # Evaluate case
    cl, cls, cms, cps, panels = solve_case(airfoil, uniform; num_panels = 80)

    @test cl       ≈  0.84188988 atol = 1e-6
    @test sum(cls) ≈  0.84073703 atol = 1e-6
    @test sum(cms) ≈ -0.26104277 atol = 1e-6
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

    @test cl       ≈  0.85736965 atol = 1e-6
    @test sum(cls) ≈  0.85976886 atol = 1e-6
    @test sum(cms) ≈ -0.29766116 atol = 1e-6
end

@testset "Trapezoidal Wing Geometry" begin
    # Define wing
    wing_right = HalfWing(chords    = [1.0, 0.6, 0.2],
                          twists    = [2.0, 0.0, -0.2],
                          spans     = [5.0, 0.5],
                          dihedrals = [5., 5.],
                          sweep_LEs = [5., 5.]);
        
    # Get wing info
    b, S, c, AR = info(wing_right)
    λ           = taper_ratio(wing_right)
    wing_mac    = mean_aerodynamic_center(wing_right)

    @test b        ≈ 5.50000000                    atol = 1e-6
    @test S        ≈ 4.19939047                    atol = 1e-6
    @test c        ≈ 0.79841269                    atol = 1e-6
    @test AR       ≈ 7.20342634                    atol = 1e-6
    @test λ        ≈ 0.20000000                    atol = 1e-6
    @test wing_mac ≈ [0.42092866, 1.33432539, 0.0] atol = 1e-6
end

@testset "Vortex Lattice Method - NACA 0012 Rectangular Wing" begin
    # Define wing
    wing = Wing(foils     = Foil.(naca4((0,0,1,2)) for i ∈ 1:2),
                chords    = [0.18, 0.16],
                twists    = [0., 0.],
                spans     = [0.5,],
                dihedrals = [5.],
                sweep_LEs = [1.14])

    # Define reference values
    ρ   = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω   = [0.0, 0.0, 0.0]
    
    # Define freestream condition
    uniform = Freestream(10.0, 2.0, 2.0, Ω)

    # Evaluate stability case
    nf_coeffs, ff_coeffs, dv_coeffs = solve_stability_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5)

    # Test values
    nf_tests = [0.001189, -0.000228, 0.152203, -0.000242, -0.003486, -8.1e-5, 0.0, 0.0, 0.0]
    ff_tests = [0.00123,  -0.000271, 0.152198]
    dv_tests = [ 0.068444 -0.000046 -0.000711  0.023607  0.000337; 
                 0.010867 -0.007536  0.129968  0.021929 -0.012086; 
                 4.402229 -0.012973 -0.070654  6.833903  0.001999; 
                 0.031877 -0.013083  0.460035  0.091216 -0.039146; 
                -0.112285 -0.004631  0.105695 -0.852395 -0.007696; 
                -0.002218 -0.002115  0.008263 -0.003817  0.001079]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nf_coeffs, nf_tests) ]

    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ff_coeffs, ff_tests) ]

    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dv_coeffs, dv_tests) ]
end;