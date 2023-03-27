using AeroFuse
using Test
using StaticArrays
using ComponentArrays

@testset "Loading Airfoil Coordinates" begin
    foil = read_foil(string(@__DIR__, "/e49.dat"))

    test_foil = Foil(Float16[ 
        1.0000000 0.0000000
        0.9966400 0.0016100
        0.9873600 0.0066800
        0.9378100 0.0375600
        0.8905500 0.0593500
        0.6676800 0.1042800
        0.6220800 0.1086800
        0.3843200 0.1143100
        0.2934500 0.1088000
        0.1385800 0.0863300
        0.0573200 0.0616200
        0.0377900 0.0525200
        0.0106500 0.0350900
        0.0032500 0.0275500
        0.0002000 0.0216700
        0.0023000 0.0183600
        0.0110000 0.0167400
        0.3336100 0.0399300
        0.3901600 0.0422400
        0.4483800 0.0437000
        0.7861900 0.0330400
        0.9856000 0.0039800
        0.9963900 0.0010800
        1.0000000 0.0000000 
    ], "EPPLER 49 AIRFOIL")

    @test coordinates(foil) == coordinates(test_foil)
end

@testset "NACA-4 Doublet-Source 2D Panel Method" begin
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

@testset "Airfoil Processing and Doublet-Source 2D Panel Method" begin
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
    wing = Wing(
        chords    = [1.0, 0.6, 0.2],
        twists    = [2.0, 0.0, -0.2],
        spans     = [5.0, 0.5],
        dihedrals = [5., 5.],
        sweeps    = [5., 10.],
    #   symmetry  = true
    );

    # Get wing info
    b = span(wing)
    S = projected_area(wing)
    c = mean_aerodynamic_chord(wing)
    AR = aspect_ratio(wing)
    λ = taper_ratio(wing)
    wing_mac = mean_aerodynamic_center(wing)

    @test b        ≈ 5.50000000                    atol = 1e-6
    @test S        ≈ 4.19939047                    atol = 1e-6
    @test c        ≈ 0.79841008                    atol = 1e-6
    @test AR       ≈ 7.20342634                    atol = 1e-6
    @test λ        ≈ 0.20000000                    atol = 1e-6
    @test wing_mac ≈ [0.4218125, 2.4305755, 0.0]   atol = 1e-6
end

@testset "Geometry - Single-Section Trapezoidal Wing" begin
    # Define wing section
    wing_sec = WingSection(aspect = 6.25, area = 4.0, taper = 0.6, symmetry = true)

    # Using Wing constructor
    wing = Wing(
        foils     = fill(naca4((0,0,1,2)), 2),
        chords    = [1.0, 0.6],
        twists    = [0.0, 0.0],
        spans     = [5.0] / 2,
        dihedrals = [11.39],
        sweeps    = [0.],
        symmetry  = true, 
    );

    @test span(wing) ≈ span(wing_sec) atol = 1e-6  # Span 
    @test projected_area(wing) ≈ projected_area(wing_sec) atol = 1e-6 # Projected area 
    @test mean_aerodynamic_chord(wing) ≈ mean_aerodynamic_chord(wing_sec) atol = 1e-6 # Mean aerodynamic chord
    @test aspect_ratio(wing) ≈ aspect_ratio(wing_sec) atol = 1e-6 # Aspect ratio
    @test taper_ratio(wing) ≈ taper_ratio(wing_sec) atol = 1e-6 # Taper ratio
    @test mean_aerodynamic_center(wing) ≈ mean_aerodynamic_center(wing_sec) atol = 1e-6 # Mean aerodynamic center
end

@testset "Geometry - 3D Panel" begin
    panel = Panel3D(
        [1.0, -1., 0.0], 
        [0.0, -1., -0.0], 
        [0.0, 0.0, -0.0], 
        [1.0, 0.0, 0.0]
    )
    
    mp = midpoint(panel)
    ε = 1.
    p = SVector(mp.x + ε, mp.y + ε, mp.z + ε)
    T = get_transformation(panel)

    @test inv(T)(T(mp) + T(p)) ≈ p atol = 1e-6
    @test panel_area(panel) ≈ 1.0 atol = 1e-6
    @test normal_vector(panel) ≈ [0,0,-2.0] atol = 1e-6
end

@testset "Freestream 3D Velocity Conversion" begin
    φ, θ = 1.0, 1.0
    fs = Freestream(alpha = φ, beta = θ)
    V_test = [ cosd(θ) * cosd(φ), -sind(φ), sind(θ) * cosd(φ) ]
    V_run = velocity(fs)

    φ_test, θ_test = cartesian_to_freestream(V_run)

    @test V_run ≈ V_test atol = 1e-6
    @test φ_test ≈ φ atol = 1e-6
    @test θ_test ≈ -θ atol = 1e-6
end

@testset "Aerodynamics - Reference Values" begin
    refs = References(
        speed    = 150.0,
        area     = 30.0,
        span     = 15.0,
        chord    = 2.0,
        density  = 1.225,
        location = [0.25, 0., 0.]
    )

    @test mach_number(refs) ≈ 0.454545 atol = 1e-6
    @test reynolds_number(refs) ≈ 2.45e7 atol = 1e-6
    @test dynamic_pressure(refs) ≈ 13781.25 atol = 1e-6
end

@testset "Vortex Lattice Method (Horseshoes, Incompressible) - NACA 0012 Tapered Wing" begin
    # Define wing
    wing = Wing(
        foils     = [ naca4((0,0,1,2)) for i ∈ 1:2 ],
        chords    = [0.18, 0.16],
        twists    = [0., 0.],
        spans     = [0.25,],
        dihedrals = [5.],
        sweeps    = [1.14],
        symmetry  = true
    )

    # Define freestream and reference values
    fs = Freestream(2.0, 2.0, [0.0, 0.0, 0.0])
    refs = References(
        speed    = 1.0, 
        area     = projected_area(wing), 
        span     = span(wing), 
        chord    = mean_aerodynamic_chord(wing), 
        density  = 1.225, 
        location = [0.25 * mean_aerodynamic_chord(wing), 0., 0.]
    )

    aircraft = ComponentArray(wing = make_horseshoes(WingMesh(wing, [20], 5, span_spacing = [Sine(1); Sine()])))

    # Evaluate stability case
    system = VortexLatticeSystem(aircraft, fs, refs, false)
    dv_data = freestream_derivatives(system)

    dcf = dv_data.wing
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,3:end]

    # Test values
    nf_tests = [0.0013533, 0.0002199, 0.1159468, 0.0009056, 0.000387, -9.45e-5]
    ff_tests = [0.0014281, 0.0001743, 0.1159415]
    dv_tests = [ 0.0797418  -0.0004935 -0.0039018  0.06373702665595818     0.0004043938437192077;
                 0.0066032   0.0024036  0.0809152  0.030210659199330615   -0.012157465112956763;
                 3.4110144  -0.0060237 -0.1856945  5.1791598961593985      0.006557296956889664;
                 0.0084356   0.0064222  0.2850743  0.11097309781645154    -0.03857231877442366;
                -0.0088809   0.0004686  0.0618588 -0.6604736234860606     -0.005118057094076678;
                -0.0009645  -0.0018438  0.0054143 -0.004469695166159752    0.001369403002390149]
    
    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end

# Include common aircraft definition
include(string(@__DIR__, "/aircraft_definition.jl"))

@testset "Vortex Lattice Method (Horseshoes, Compressible) - Vanilla Aircraft" begin
    aircraft = ComponentArray(
        wing  = make_horseshoes(wing_mesh),
        htail = make_horseshoes(htail_mesh),
        vtail = make_horseshoes(vtail_mesh),
    )
    ## Stability case
    system = VortexLatticeSystem(aircraft, fs, refs, true)
    dv_data = freestream_derivatives(system)

    dcf = dv_data.aircraft
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,2:end]

    nf_tests = [0.0004394, -0.0096616, 0.0710165, -0.0027487, 0.0643144, 0.0063322]
    ff_tests = [0.0005934, -0.0097369, 0.0709867]
    dv_tests = [ 0.000303126   0.0318917    0.00657018   0.00105799    0.0840928  -0.000608665
                -0.00222327    0.00632911  -0.542542     0.31718      -0.0952224  -1.07032
                0.0285508     5.164033      0.0671273   -0.0606896    14.2909359      0.0398133
                -0.000318489   0.0468793   -0.15929      0.514803      0.184925   -0.11513
                0.0205454    -1.5507      -0.0366944   -0.0862627   -29.7754714     -0.149386
                0.00161937    0.0117921    0.35447      0.0190155     0.119573    0.825448 ]

    # Nearfield coefficients test
    [ @test nf_c ≈ nf_t atol = 1e-6 for (nf_c, nf_t) in zip(nfs, nf_tests) ]
    # Farfield coefficients test
    [ @test ff_c ≈ ff_t atol = 1e-6 for (ff_c, ff_t) in zip(ffs, ff_tests) ]
    # # Stability derivatives' coefficients test
    [ @test dv_c ≈ dv_t atol = 1e-6 for (dv_c, dv_t) in zip(dvs, dv_tests) ]
end

@testset "Vortex Lattice Method (Rings, Compressible) - Vanilla Aircraft" begin
    aircraft = ComponentArray(
        wing  = make_vortex_rings(wing_mesh),
        htail = make_vortex_rings(htail_mesh),
        vtail = make_vortex_rings(vtail_mesh),
    )

    ## Stability case
    dv_data = freestream_derivatives(aircraft, fs, refs, name = :aircraft, compressible = true)

    dcf = dv_data.aircraft
    nfs = @views dcf[1:6,1]
    ffs = @views dcf[7:9,1]
    dvs = @views dcf[1:6,2:end]

    nf_tests = [0.0004485, -0.0102081, 0.0706159, -0.0028898, 0.0653328, 0.0067422]
    ff_tests = [0.0006049, -0.0102852, 0.0705857]
    dv_tests = [ 0.000306494   0.0316461    0.00767591   0.00118087   0.0837572   0.000839109
                -0.00232925    0.00389612  -0.565213     0.315484    -0.0809676  -1.13528256
                 0.0283188     5.1464928    0.0345146   -0.04648     14.2846332   0.0138955
                -0.000401224   0.0417969   -0.167955     0.515156     0.16723    -0.120098
                 0.021139     -1.502605     0.0895636   -0.0791555  -29.800452   -0.0345743
                 0.00166363    0.0116769    0.370111     0.0211191    0.102458    0.882128 ]

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