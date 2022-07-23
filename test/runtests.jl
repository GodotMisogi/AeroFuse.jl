using AeroMDAO
using Test

@testset "NACA-4 Doublet-Source Panel Method" begin
    # Define airfoil
    airfoil = (naca4)((0,0,1,2))

    # Define uniform flow
    uniform = Uniform2D(1., 0.)

    # Evaluate case
    sys_1 = solve_case(airfoil, uniform; num_panels = 80)
    cl_1  = lift_coefficient(sys_1)
    cls_1 = surface_coefficients(sys_1)[1]

    # α = 5ᵒ
    uniform = Uniform2D(1., 5.)
    sys_2 = solve_case(airfoil, uniform; num_panels = 80)
    cl_2  = lift_coefficient(sys_2)
    cls_2 = surface_coefficients(sys_2)[1]

    @test cl_1       ≈ 0.0       atol = 1e-6
    @test cl_2       ≈ 0.5996184 atol = 1e-6
    @test sum(cls_2) ≈ 0.6007449 atol = 1e-6
end

@testset "Airfoil Processing and Doublet-Source Panel Method" begin
    # Import and read airfoil coordinates
    coo_foil = naca4((2,4,1,2))

    # Cosine spacing
    cos_foil = cosine_interpolation(coo_foil, 61)

    # Split airfoil
    up, low  = split_surface(cos_foil)

    # Convert coordinates to Kulfan CST variables
    num_dv   = 4
    alpha_u  = coordinates_to_CST(up, num_dv)
    alpha_l  = coordinates_to_CST(low, num_dv)

    # Generate same airfoil using Kulfan CST parametrisation
    cst_foil = kulfan_CST(alpha_u, alpha_l, (0., 0.), (0., 0.))

    # Test coefficients
    uniform  = Uniform2D(1., 5.)

    sys_coo            = solve_case(coo_foil, uniform; num_panels = 80)
    cl_coo             = lift_coefficient(sys_coo)
    cls_coo, cms_coo   = surface_coefficients(sys_coo)[1:2]
    @test cl_coo       ≈  0.83220516 atol = 1e-6
    @test sum(cls_coo) ≈  0.83291636 atol = 1e-6
    @test sum(cms_coo) ≈ -0.25899389 atol = 1e-6

    sys_cos            = solve_case(cos_foil, uniform; num_panels = 80)
    cl_cos             = lift_coefficient(sys_cos)
    cls_cos, cms_cos   = surface_coefficients(sys_cos)[1:2]
    @test cl_cos       ≈  0.83178821 atol = 1e-6
    @test sum(cls_cos) ≈  0.83269773 atol = 1e-6
    @test sum(cms_cos) ≈ -0.25889408 atol = 1e-6

    sys_cst            = solve_case(cst_foil, uniform; num_panels = 80)
    cl_cst             = lift_coefficient(sys_cst)
    cls_cst, cms_cst   = surface_coefficients(sys_cst)[1:2]
    @test cl_cst       ≈  0.83381613 atol = 1e-6
    @test sum(cls_cst) ≈  0.83408259 atol = 1e-6
    @test sum(cms_cst) ≈ -0.25986701 atol = 1e-6
end

@testset "Geometry - Two-Section Trapezoidal Wing" begin
    # Define wing
    wing_right = Wing(chords    = [1.0, 0.6, 0.2],
                          twists    = [2.0, 0.0, -0.2],
                          spans     = [5.0, 0.5],
                          dihedrals = [5., 5.],
                          sweeps    = [5., 5.],
                        #   symmetry  = true
                        );

    # Get wing info
    b        = span(wing_right)
    S        = projected_area(wing_right)
    c        = mean_aerodynamic_chord(wing_right)
    AR       = aspect_ratio(wing_right)
    λ        = taper_ratio(wing_right)
    wing_mac = mean_aerodynamic_center(wing_right)

    @test b        ≈ 5.50000000                    atol = 1e-6
    @test S        ≈ 4.19939047                    atol = 1e-6
    @test c        ≈ 0.79841008                    atol = 1e-6
    @test AR       ≈ 7.20342634                    atol = 1e-6
    @test λ        ≈ 0.20000000                    atol = 1e-6
    @test wing_mac ≈ [0.4209310, 1.3343524, 0.0] atol = 1e-6
end

@testset "Vortex Lattice Method (Incompressible) - NACA 0012 Tapered Wing" begin
    # Define wing
    wing = Wing(foils     = [ naca4((0,0,1,2)) for i ∈ 1:2 ],
                chords    = [0.18, 0.16],
                twists    = [0., 0.],
                spans     = [0.25,],
                dihedrals = [5.],
                sweeps    = [1.14],
                symmetry  = true
            )

    # Define freestream and reference values
    fs   = Freestream(2.0, 2.0, [0.0, 0.0, 0.0])
    refs = References(speed    = 1.0, 
                      area     = projected_area(wing), 
                      span     = span(wing), 
                      chord    = mean_aerodynamic_chord(wing), 
                      density  = 1.225, 
                      location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.])

    aircraft = ComponentArray(wing = make_horseshoes(WingMesh(wing, [20], 5, span_spacing = [Sine(1); Sine()])))

    # Evaluate stability case
    dv_data = solve_case_derivatives(aircraft, fs, refs)

    dcf = dv_data.wing
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,3:end]

    # Test values
    nf_tests = [0.0013533, 0.0002199, 0.1159468, 0.0009056, 0.000387, -9.45e-5]
    ff_tests = [0.0014281, 0.0001743, 0.1159415]
    dv_tests = [ 0.0797419  -0.0004936  -0.0039018   0.0637380   0.0004044
                 0.0066033   0.0024036   0.0809152   0.0302107  -0.0121575
                 3.4110222  -0.0060238  -0.185695    5.179196    0.0065573
                 0.0084357   0.0064223   0.285074    0.1109742  -0.0385723
                -0.0088809   0.0004686   0.0618588  -0.6604790  -0.0051181
                -0.0009645  -0.0018438   0.0054144  -0.0044697   0.0013694]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end

@testset "Vortex Lattice Method (Incompressible) - Vanilla Aircraft" begin
    ## Wing
    wing = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [1.0, 0.6],
        twists    = [0.0, 0.0],
        spans     = [5.0] / 2,
        dihedrals = [11.39],
        sweeps    = [0.],
        symmetry  = true, 
    );

    # Horizontal tail
    htail = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.25] / 2,
        dihedrals = [0.],
        sweeps    = [6.39],
        position  = [4., 0, 0],
        angle     = -2.,
        axis      = [0., 1., 0.],
        symmetry  = true
    )

    # Vertical tail
    vtail = Wing(
        foils     = fill(naca4((0,0,0,9)), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.0],
        dihedrals = [0.],
        sweeps    = [7.97],
        position  = [4., 0, 0],
        angle     = 90.,
        axis      = [1., 0., 0.],
    )

    ## Assembly
    wing_mesh = WingMesh(wing, [32], 10; span_spacing = Cosine())
    htail_mesh = WingMesh(htail, [12],  6; span_spacing = Cosine())
    vtail_mesh = WingMesh(vtail, [5],  6; span_spacing = Cosine())

    aircraft = ComponentArray(
        wing  = make_horseshoes(wing_mesh),
        htail = make_horseshoes(htail_mesh),
        vtail = make_horseshoes(vtail_mesh),
    )

    ## Reference quantities
    fs = Freestream(
        alpha    = 1.0, 
        beta     = 1.0, 
        omega    = zeros(3)
    )
                         
    refs = References(
        speed    = 1.0,
        area     = projected_area(wing),
        span     = span(wing),
        chord    = mean_aerodynamic_chord(wing),
        density  = 1.225,
        location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    )

    ## Stability case
    dv_data = solve_case_derivatives(aircraft, fs, refs);

    dcf = dv_data.aircraft
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,3:end]

    nf_tests = [0.0003809, -0.0092002, 0.0653607, -0.0026762, 0.0601839, 0.0059991]
    ff_tests = [0.0005143, -0.0092706, 0.0653322]
    dv_tests = [ 0.0268317   0.0059901   0.0012994    0.0767768  -0.0023994
                 0.0048546  -0.517466    0.30401     -0.0744164  -0.899784
                 4.7952172   0.0598145  -0.0557912   11.9681911   0.0348196
                 0.0413913  -0.155083    0.489466     0.1469693  -0.0990389
                -1.5608632  -0.0360026  -0.0616955  -25.5170944  -0.124964
                 0.010348    0.336414    0.0155232    0.0906727   0.693076 ]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end

@testset "Vortex Lattice Method (Compressible) - Vanilla Aircraft" begin
    ## Wing
    wing = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [1.0, 0.6],
        twists    = [0.0, 0.0],
        spans     = [5.0] / 2,
        dihedrals = [11.39],
        sweeps    = [0.],
        symmetry  = true, 
    );

    # Horizontal tail
    htail = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.25] / 2,
        dihedrals = [0.],
        sweeps    = [6.39],
        position  = [4., 0, 0],
        angle     = -2.,
        axis      = [0., 1., 0.],
        symmetry  = true
    )

    # Vertical tail
    vtail = Wing(
        foils     = fill(naca4((0,0,0,9)), 2),
        chords    = [0.7, 0.42],
        twists    = [0.0, 0.0],
        spans     = [1.0],
        dihedrals = [0.],
        sweeps    = [7.97],
        position  = [4., 0, 0],
        angle     = 90.,
        axis      = [1., 0., 0.],
    )

    ## Assembly
    wing_mesh = WingMesh(wing, [32], 10; span_spacing = Cosine())
    htail_mesh = WingMesh(htail, [12], 6; span_spacing = Cosine())
    vtail_mesh = WingMesh(vtail, [5], 6; span_spacing = Cosine())

    aircraft = ComponentArray(
        wing  = make_horseshoes(wing_mesh),
        htail = make_horseshoes(htail_mesh),
        vtail = make_horseshoes(vtail_mesh),
    )

    ## Reference quantities
    fs      = Freestream(alpha    = 1.0, 
                         beta     = 1.0, 
                         omega    = zeros(3))
                         
    refs    = References(speed    = 150.0,
                         area     = projected_area(wing),
                         span     = span(wing),
                         chord    = mean_aerodynamic_chord(wing),
                         density  = 1.225,
                         location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.])

    ## Stability case
    dv_data = solve_case_derivatives(aircraft, fs, refs);

    dcf = dv_data.aircraft
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,3:end]

    nf_tests = [0.0004394, -0.0096616, 0.0710165, -0.0027487, 0.0643144, 0.0063322]
    ff_tests = [0.0005934, -0.0097369, 0.0709867]
    dv_tests = [ 0.0318917  0.0065702  -0.6380811    -1.3022877    1.8906515;
                 0.0063291 -0.542542   47.3510327   -16.8276434 -159.4190293;
                 5.1640333  0.0671273  -9.0759321   2150.559606    5.2457761;
                 0.0468793 -0.1592905  77.192654     27.2895215  -15.6091835;
                -1.5506993 -0.0366944 -13.4122289 -4464.6998256  -20.6719332;
                 0.0117921  0.3544702   1.1794656    19.5810734  123.8328465]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # # Stability derivatives' coefficients test
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