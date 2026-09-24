# %%
using AeroFuse
using Test
using StaticArrays
using LinearAlgebra
import AeroFuse.VortexLattice: AbstractVortex, VortexLatticeSystem
import AeroFuse.VortexLattice: bound_leg_velocity

struct NonWakeTestElement <: AbstractPotentialFlowElement end

@testset "Potential-flow API compatibility" begin
    @test AeroFuse.VortexLattice === AeroFuse.PotentialFlow
    @test VortexLatticeSystem === PotentialFlowSystem
    @test AbstractVortex === AbstractPotentialFlowElement
    @test bound_leg_velocity === AeroFuse.PotentialFlow.bound_leg_velocity
    @test !AeroFuse.PotentialFlow.has_wake(NonWakeTestElement())
end

# %%
@testset "Coupled potential-flow post-processing" begin
    wing = Wing(
        foils = fill(naca4(0, 0, 1, 2), 2),
        chords = [0.6, 0.4], spans = [1.0], symmetry = true,
    )
    mesh = WingMesh(wing, [4], 2)
    vortices = elements(mesh, Horseshoe())
    fuselage = HyperEllipseFuselage(
        radius = 0.15, length = 2.0, c_nose = 2, c_rear = 2,
        position = [-0.5, 0.0, -0.5],
    )
    panels = elements(fuselage, SourcePanel(); n_secs = 4, n_circ = 7)
    lines = make_fuselage_line(fuselage; n = 4)
    fs = Freestream(alpha = 5.0, beta = 2.0)
    refs = References(speed = 10.0, location = [0.1, 0.0, 0.0])

    for (name, nonlifting) in ((:skin, panels), (:slender, lines))
        aircraft = ComponentVector(; wing = vortices, name => nonlifting)
        system = PotentialFlowSystem(aircraft, fs, refs)
        @test system.vortices === system.elements
        @test system.circulations === system.strengths
        @test all(name -> hasproperty(system, name),
                  (:elements, :strengths, :vortices, :circulations))
        @test fieldnames(typeof(system))[1:2] == (:elements, :strengths)
        @test all(isfinite, system.strengths)
        @test system.influence_matrix * (system.strengths / refs.speed) ≈
            system.boundary_vector
        @test VortexLatticeSystem(aircraft, fs, refs).strengths ≈ system.strengths
        @test farfield_forces(system)[name] == zeros(SVector{3,Float64})
        @test AeroFuse.PotentialFlow.has_wake(first(vortices))

        expected = name === :skin ? body_forces(system, name) :
            AeroFuse.PotentialFlow.fuselage_munk_forces(system, name)
        @test any(force -> norm(force) > 0, expected)
        @test surface_forces(system)[name] ≈ expected
        for axes in (Geometry(), Body(), Wind(), Stability())
            forces, moments = surface_dynamics(system; axes)
            @test surface_forces(system; axes) ≈ forces
            @test surface_moments(system; axes) ≈ moments
        end

        renamed = PotentialFlowSystem(
            ComponentVector(wing = vortices, renamed = nonlifting), fs, refs)
        @test renamed.strengths.renamed ≈ system.strengths[name]
        @test surface_forces(renamed).renamed ≈ expected
    end

    @test_throws ArgumentError PotentialFlowSystem(
        ComponentVector(mixed = AbstractPotentialFlowElement[
            first(vortices), first(panels)]), fs, refs)
    @test_throws ArgumentError PotentialFlowSystem(
        ComponentVector(fuse = lines[1:1]), fs, refs)
    @test_throws ArgumentError PotentialFlowSystem(
        ComponentVector(fuse = [first(lines), first(lines)]), fs, refs)
    @test_throws ArgumentError PotentialFlowSystem(
        ComponentVector(fuse = reverse(lines)), fs, refs)
    @test_throws ArgumentError PotentialFlowSystem(
        ComponentVector(wing = [1.0]), fs, refs)
    @test_throws ArgumentError PotentialFlowSystem(vortices, fs, refs)

    A = zeros(length(vortices), length(vortices))
    @test AeroFuse.PotentialFlow.influence_matrix!(A, vec(vortices)) === A
    @test A ≈ AeroFuse.PotentialFlow.influence_matrix(vec(vortices))
    @test_throws DimensionMismatch AeroFuse.PotentialFlow.influence_matrix!(
        zeros(1, 1), vec(vortices))

    rings = elements(mesh, VortexRing())
    @test AeroFuse.PotentialFlow.has_wake(first(rings))
    @test all(isfinite, PotentialFlowSystem(
        ComponentVector(wing = rings), fs, refs, true).strengths)

    prop = PropellerDisk(
        center = [-1.0, 0.0, 0.0], axis = [1.0, 0.0, 0.0],
        radius = 1.5, thrust = 10.0, turn = 0.2,
    )
    plain = PotentialFlowSystem(ComponentVector(wing = vortices), fs, refs)
    blown = PotentialFlowSystem(ComponentVector(wing = vortices), fs, refs;
                                slipstream = prop)
    @test all(isfinite, blown.strengths)
    @test !isapprox(blown.strengths, plain.strengths)
    @test surface_forces(blown) ≈ first(surface_dynamics(blown))
end

# %%
using KernelAbstractions: CPU

@testset "Potential-flow compute backends" begin
    wing = Wing(
        foils = fill(naca4(0, 0, 1, 2), 2),
        chords = [0.6, 0.4], spans = [1.0], dihedrals = [5.0], symmetry = true,
    )
    mesh = WingMesh(wing, [8], 4)
    fuselage = HyperEllipseFuselage(
        radius = 0.15, length = 2.0, c_nose = 2, c_rear = 2,
        position = [-0.5, 0.0, -0.5],
    )
    fs = Freestream(alpha = 5.0, beta = 2.0, omega = [0.1, 0.2, 0.0])
    refs = References(speed = 10.0, location = [0.1, 0.0, 0.0])
    prop = PropellerDisk(
        center = [-1.0, 0.0, 0.0], axis = [1.0, 0.0, 0.0],
        radius = 1.5, thrust = 10.0, turn = 0.2,
    )

    @test Base.get_extension(AeroFuse, :AeroFuseKernelAbstractionsExt) !== nothing

    cases = (
        skin    = (ComponentVector(wing = elements(mesh, Horseshoe()), body = elements(fuselage, SourcePanel(); n_secs = 4, n_circ = 7)), nothing),
        slender = (ComponentVector(wing = elements(mesh, Horseshoe()), fuse = make_fuselage_line(fuselage; n = 4)), nothing),
        rings   = (ComponentVector(wing = elements(mesh, VortexRing())), nothing),
        blown   = (ComponentVector(wing = elements(mesh, Horseshoe())), prop),
    )

    for (name, (aircraft, slipstream)) in pairs(cases), compressible in (false, true)
        ref = PotentialFlowSystem(aircraft, fs, refs, compressible; slipstream)
        sys = PotentialFlowSystem(aircraft, fs, refs, compressible; slipstream, backend = CPU())

        @test sys.backend === CPU()
        @test keys(sys.strengths) == keys(ref.strengths)
        @test sys.influence_matrix ≈ ref.influence_matrix rtol = 1e-12
        @test sys.boundary_vector ≈ ref.boundary_vector rtol = 1e-12
        @test sys.strengths ≈ ref.strengths rtol = 1e-10
        @test reduce(vcat, surface_velocities(sys)) ≈ reduce(vcat, surface_velocities(ref)) rtol = 1e-10
        @test sum(surface_forces(sys)) ≈ sum(surface_forces(ref)) rtol = 1e-10
        @test sum(surface_moments(sys)) ≈ sum(surface_moments(ref)) rtol = 1e-10
        @test farfield(sys) ≈ farfield(ref) rtol = 1e-10
        name === :skin && @test reduce(vcat, body_surface_velocities(sys)) ≈ reduce(vcat, body_surface_velocities(ref)) rtol = 1e-10
    end

    # Float32 kernels (used by GPU backends without Float64) must stay finite and accurate:
    # guards the cancellation-free vortex kernel forms.
    ext = Base.get_extension(AeroFuse, :AeroFuseKernelAbstractionsExt)
    aircraft = first(cases.skin)
    wind = map(el -> AeroFuse.PotentialFlow.geometry_to_wind_axes(el, fs), vec(collect(aircraft)))
    A64 = AeroFuse.PotentialFlow.influence_matrix(wind)
    A32 = AeroFuse.PotentialFlow.influence_matrix([ ext.convert_element(Float32, el) for el in wind ])
    @test eltype(A32) == Float32
    @test all(isfinite, A32)
    @test norm(A32 - A64) / norm(A64) < 1e-5

    # Points on a trailing line's axis (upstream: |r| + r·u = 0 exactly) induce no velocity,
    # rather than NaN from the core term 0/0.
    u = SVector(1f0, 0f0, 0f0)
    for r in (SVector(-1f0, 0f0, 0f0), SVector(-1.0, 0.0, 0.0))
        @test AeroFuse.PotentialFlow.trailing_leg_velocity(r, one(eltype(r)), u, zero(eltype(r))) == zero(r)
        @test AeroFuse.PotentialFlow.trailing_leg_velocity(r, one(eltype(r)), u) == zero(r)
    end
end
