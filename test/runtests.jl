using AeroMDAO
using Test
using ComponentArrays

@testset "NACA-4 Doublet-Source Panel Method" begin
    # Define airfoil
    airfoil = (Foil ∘ naca4)((0,0,1,2))

    # Define uniform flow
    uniform = Uniform2D(1., 0.)

    # Evaluate case
    cl_1, cls_1, cms_1, cps_1, panels = solve_case(airfoil, uniform; num_panels = 80)

    # α = 5ᵒ
    uniform = Uniform2D(1., 5.)
    cl_2, cls_2, cms_2, cps_2, panels = solve_case(airfoil, uniform; num_panels = 80)

    @test cl_1       ≈ 0.0       atol = 1e-6
    @test cl_2       ≈ 0.5996184 atol = 1e-6
    @test sum(cls_2) ≈ 0.6007449 atol = 1e-6
end

@testset "Airfoil Processing and Doublet-Source Panel Method" begin
    # Import and read airfoil coordinates
    coo_foil  = naca4(2,4,1,2)

    # Cosine spacing
    cos_foil = cosine_foil(coo_foil, 61)

    # Split airfoil
    up, low  = split_foil(cos_foil)

    # Convert coordinates to Kulfan CST variables
    num_dv   = 4
    alpha_u = coords_to_CST(up, num_dv)
    alpha_l = coords_to_CST(low, num_dv)

    # Generate same airfoil using Kulfan CST parametrisation
    cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), (0., 0.))

    # Test coefficients
    uniform  = Uniform2D(1., 5.)

    cl_coo, cls_coo, cms_coo = solve_case(Foil(coo_foil), uniform; num_panels = 80)[1:3]
    @test cl_coo       ≈  0.83220516 atol = 1e-6
    @test sum(cls_coo) ≈  0.83291636 atol = 1e-6
    @test sum(cms_coo) ≈ -0.25899389 atol = 1e-6

    cl_cos, cls_cos, cms_cos = solve_case(Foil(cos_foil), uniform; num_panels = 80)[1:3]
    @test cl_cos       ≈  0.83178821 atol = 1e-6
    @test sum(cls_cos) ≈  0.83269773 atol = 1e-6
    @test sum(cms_cos) ≈ -0.25889408 atol = 1e-6

    cl_cst, cls_cst, cms_cst = solve_case(Foil(cst_foil), uniform; num_panels = 80)[1:3]
    @test cl_cst       ≈  0.83381613 atol = 1e-6
    @test sum(cls_cst) ≈  0.83408259 atol = 1e-6
    @test sum(cms_cst) ≈ -0.25986701 atol = 1e-6
end

@testset "Geometry - Two-Section Trapezoidal Wing" begin
    # Define wing
    wing_right = HalfWing(chords    = [1.0, 0.6, 0.2],
                          twists    = [2.0, 0.0, -0.2],
                          spans     = [5.0, 0.5],
                          dihedrals = [5., 5.],
                          LE_sweeps = [5., 5.]);

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
                LE_sweeps = [1.14])

    # Define reference values
    ρ   = 1.225
    ref = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    Ω   = [0.0, 0.0, 0.0]

    # Define freestream condition
    uniform = Freestream(10.0, 2.0, 2.0, Ω)

    # Evaluate stability case
    nf_coeffs, ff_coeffs, dv_coeffs = solve_stability_case(wing, uniform; rho_ref = ρ, r_ref = ref, span_num = 20, chord_num = 5)

    # Test values
    nf_tests = [0.001189, -0.000228, 0.152203, -0.000242, -0.003486, -8.1e-5]
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
end

@testset "Vortex Lattice Method - Vanilla Aircraft" begin
    ## Wing
    wing = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                chords    = [1.0, 0.6],
                twists    = [0.0, 0.0],
                spans     = [5.0],
                dihedrals = [11.39],
                LE_sweeps = [0.]);

    # Horizontal tail
    htail = Wing(foils     = Foil.(fill(naca4((0,0,1,2)), 2)),
                 chords    = [0.7, 0.42],
                 twists    = [0.0, 0.0],
                 spans     = [1.25],
                 dihedrals = [0.],
                 LE_sweeps = [6.39],
                 position  = [4., 0, 0],
                 angle     = -2.,
                 axis      = [0., 1., 0.])

    # Vertical tail
    vtail = HalfWing(foils     = Foil.(fill(naca4((0,0,0,9)), 2)),
                     chords    = [0.7, 0.42],
                     twists    = [0.0, 0.0],
                     spans     = [1.0],
                     dihedrals = [0.],
                     LE_sweeps = [7.97],
                     position  = [4., 0, 0],
                     angle     = 90.,
                     axis      = [1., 0., 0.])

    ## Assembly
    wing_panels , wing_normals  = panel_wing(wing, 16, 10; spacing = "cosine")
    htail_panels, htail_normals = panel_wing(htail, 6,  6; spacing = "cosine")
    vtail_panels, vtail_normals = panel_wing(vtail, 5,  6; spacing = "cosine")

    aircraft = ComponentArray(
                              wing  = Horseshoe.(wing_panels , wing_normals),
                              htail = Horseshoe.(htail_panels, htail_normals),
                              vtail = Horseshoe.(vtail_panels, vtail_normals)
                             )

    ## Reference quantities
    ac_name = :aircraft
    S, b, c = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
    ρ       = 1.225
    ref     = [0.25c, 0., 0.]
    V, α, β = 1.0, 1.0, 1.0
    Ω       = [0.0, 0.0, 0.0]
    fs      = Freestream(V, α, β, Ω)
    refs    = References(S, b, c, ρ, ref)

    ## Stability case
    dv_data = solve_stability_case(aircraft, fs, refs;
                                   name = ac_name);

    dcf = dv_data[ac_name]
    nfs = dcf.NF
    ffs = dcf.FF
    dvs = dcf.dNF

    nf_tests = [0.000258, -0.006642, 0.074301, -0.003435, 0.075511, 0.001563]
    ff_tests = [0.000375, -0.006685, 0.074281]
    dv_tests = [ 0.016795  0.003460  0.003761   0.093303 -0.000674;
                -0.000863 -0.374410  0.403476   0.000630 -0.253848;
                 5.749765  0.046649 -0.01346   15.571205  0.020396;
                 0.022674 -0.196605  0.660392   0.099065 -0.039688;
                -2.70367  -0.132928  0.070111 -37.372278 -0.064439;
                 0.002034  0.087382  0.014991   0.005840  0.091088]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end

@testset "Structures - Euler-Bernoulli Beam Elastic Stiffness" begin
    # Deflection stiffness matrix
    K = bending_stiffness_matrix([1., 1.], [1., 1.], [2., 2.], :z)

    ## 1. Fixed hinged beam subjected to force and moment at the center
    A = K[[3,4,6],[3,4,6]]  # v2, φ2, φ3
    b = [-1000, 1000, 0]    # F2, M2, M3

    x = A \ b

    ## Forces
    F1 = K * [ 0.; 0.; x[1:2]; 0.; x[3] ]

    ## 2. Propped cantilever beam with force at one end
    A = K[[1,2,4],[1,2,4]] # v1, φ1, φ2
    b = [10, 0, 0]

    x = A \ b

    ## Forces
    F2 = K * [ x[1:2]; 0.; x[3]; 0.; 0. ]

    @test F1 ≈ [968.75, 875., -1e3, 1e3, 31.25, 0.] atol = 1e-6
    @test F2 ≈ [10., 0., -25., 0., 15, -10.] atol = 1e-6
end

@testset "Structures - Euler-Bernoulli Beam Axial Stiffness" begin
    # Axial stiffness matrix
    J = axial_stiffness_matrix([1., 1., 1.], [1., 1., 1.], [2., 2., 2.])

    ## 1. ???
    A = J[[1,2],[1,2]] # ψ1, ψ2
    b = [-1000, 1000]  # R2, R2

    x = A \ b

    M = J * [ x; zeros(2) ]

    @test M ≈ [-1000., 1000., 0., 0.] atol = 1e-6
end;