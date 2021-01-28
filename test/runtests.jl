using AeroMDAO
using Test

@testset "Kulfan CST Doublet-Source Panel Method" begin
    alpha_u = [0.2, 0.3, 0.2, 0.15, 0.2]
    alpha_l = [-0.2, -0.1, -0.1, -0.001]
    dzs     = (1e-4, 1e-4)

    airfoil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, dzs, 0.0, 60);
    
    uniform = Uniform2D(1., 5.)
    cl = solve_case(airfoil, uniform, 100)

    @test isapprox(cl, 0.8911785, atol=1e-6)
end

@testset "Airfoil Processing" begin
    foilpath = joinpath((dirname ∘ dirname ∘ pathof)(AeroMDAO), "test/CRM.dat")
    coords = read_foil(foilpath)

    cos_foil = cosine_foil(coords, 51)

    up, low = split_foil(cos_foil)

    num_dv = 4
    alpha_u, alpha_l = coords_to_CST(up, num_dv), coords_to_CST(low, num_dv)

    cst_foil = (Foil ∘ kulfan_CST)(alpha_u, alpha_l, (1e-4, -1e-4), 0.0)

    uniform = Uniform2D(1., 5.)
    cl = solve_case(cst_foil, uniform, 100)

    @test isapprox(cl, 0.2486895, atol=1e-5)
end

@testset "NACA-4 Vortex Lattice Method" begin
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
    uniform = Freestream(10.0, 5.0, 5.0, Ω)
    nf_coeffs, ff_coeffs, horseshoe_panels, camber_panels, horseshoes, Γs = solve_case(wing, uniform, ref, span_num = 5, chord_num = 5) 

    CL_nf, CDi_nf, CY_nf, Cl_nf, Cm_nf, Cn_nf, p_b_nf, q_b_nf, r_b_nf = nf_coeffs
    CL_ff, CDi_ff, CY_ff, Cl_ff, Cm_ff, Cn_ff, p_b_ff, q_b_ff, r_b_ff = ff_coeffs
    
    @test isapprox(CL_ff, 0.6725620, atol=1e-5)
    @test isapprox(CDi_ff, 0.0011807, atol=1e-5)
    @test isapprox(CY_ff, -0.0037512, atol=1e-5)
    @test isapprox(Cl_nf, 0.0073144, atol=1e-5)
    @test isapprox(Cm_nf, -0.1404512, atol=1e-5)
    @test isapprox(Cn_nf, 0.0014325, atol=1e-5)
end;