# Metal (Apple GPU) backend tests for the potential-flow solver.
# Run from the repository root:
#   julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
#   julia --project=test/gpu test/gpu/runtests.jl

# %%
using AeroFuse
using ComponentArrays
using KernelAbstractions
using LinearAlgebra
using Metal
using Test

include(joinpath(@__DIR__, "..", "aircraft_definition.jl"))

fuselage = HyperEllipseFuselage(
    radius   = 0.3,
    length   = 6.,
    x_a      = 0.2,
    x_b      = 0.7,
    c_nose   = 1.6,
    c_rear   = 1.3,
    d_nose   = -0.2,
    d_rear   = 0.2,
    position = [-1., 0., 0.],
)

relerr(a, b) = norm(a - b) / norm(b)
flatten(vs) = reduce(vcat, vs)

# Float32 on Metal. Assembly is elementwise, so the AIC matches to ~ε₃₂. Solve-derived quantities
# carry κ(AIC)·ε₃₂ error: κ ≈ 7e1 for lifting surfaces alone but ≈ 6e4 with coarse source-panel
# bodies, giving up to ~3e-4 relative error on strengths and coefficients.
const RTOL_AIC   = 1e-5
const RTOL_SOLVE = 1e-3

# %%
@testset "Metal backend - $name" for (name, aircraft) in [
        "horseshoes + rings + source panels" => ComponentVector(
            wing  = elements(wing_mesh, Horseshoe()),
            htail = elements(htail_mesh, VortexRing()),
            vtail = elements(vtail_mesh, Horseshoe()),
            body  = elements(fuselage, SourcePanel(); n_secs = 6, n_circ = 8),
        ),
        "horseshoes + fuselage line" => ComponentVector(
            wing = elements(wing_mesh, Horseshoe()),
            fuse = make_fuselage_line(fuselage; n = 6),
        ),
    ]
    @test Metal.functional()

    for compressible in (false, true)
        ref = PotentialFlowSystem(aircraft, fs, refs, compressible)
        sys = PotentialFlowSystem(aircraft, fs, refs, compressible; backend = MetalBackend())

        @test sys.influence_matrix isa MtlMatrix{Float32}
        @test eltype(sys.strengths) == Float64
        @test keys(sys.strengths) == keys(ref.strengths)

        @test relerr(Array(sys.influence_matrix), ref.influence_matrix) < RTOL_AIC
        @test relerr(sys.strengths, ref.strengths) < RTOL_SOLVE
        @test relerr(flatten(surface_velocities(sys)), flatten(surface_velocities(ref))) < RTOL_SOLVE

        CFs, CMs         = surface_coefficients(sys)
        CFs_ref, CMs_ref = surface_coefficients(ref)
        @test relerr(sum(CFs), sum(CFs_ref)) < RTOL_SOLVE
        @test relerr(sum(CMs), sum(CMs_ref)) < RTOL_SOLVE
        @test relerr(farfield(sys), farfield(ref)) < RTOL_SOLVE
    end
end

@testset "Metal backend - body surface velocities" begin
    aircraft = ComponentVector(
        wing = elements(wing_mesh, Horseshoe()),
        body = elements(fuselage, SourcePanel(); n_secs = 6, n_circ = 8),
    )
    ref = PotentialFlowSystem(aircraft, fs, refs)
    sys = solve_case(aircraft, fs, refs; backend = MetalBackend())

    @test relerr(flatten(body_surface_velocities(sys)), flatten(body_surface_velocities(ref))) < RTOL_SOLVE
end
