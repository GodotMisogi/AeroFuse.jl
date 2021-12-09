using Revise
using AeroMDAO
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using ForwardDiff, ReverseDiff, Zygote
using ProtoStructs

## Foil case
#===================================================================#

x_coords = [ 1.0  0.0
             0.5  0.5
             0.0  0.0
             0.5 -0.5
             1.0  0.0 ]

## ForwardDiff (PASSES)
ForwardDiff.gradient(arc_length ∘ Foil, x_coords)

## ReverseDiff (PASSES)
ReverseDiff.gradient(arc_length ∘ Foil, x_coords)

## Zygote (PASSES)
Zygote.gradient(arc_length ∘ Foil, x_coords)


## HalfWingTest
#===================================================================#

# @proto 
struct HalfWingTest{T <: Real}
    chords    :: Vector{T}
    twists    :: Vector{T}
    spans     :: Vector{T}
    dihedrals :: Vector{T}
    sweeps    :: Vector{T}
end

HalfWingTest(chords :: AbstractVector{T}, twists :: AbstractVector{T}, spans :: AbstractVector{T}, dihedrals :: AbstractVector{T}, sweeps :: AbstractVector{T}) where {T <: Real} = HalfWingTest{T}(chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps))

##
AeroMDAO.AircraftGeometry.mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
section_macs(wing :: HalfWingTest) = @views mean_aerodynamic_chord.(wing.chords[1:end-1], fwddiv(wing.chords))
section_projected_areas(wing :: HalfWingTest) = wing.spans .* fwdsum(wing.chords) / 2

function AeroMDAO.AircraftGeometry.mean_aerodynamic_chord(wing :: HalfWingTest)
    areas = section_projected_areas(wing)
    macs  = section_macs(wing)
    sum(macs .* areas) / sum(areas)
end

wingy(x, n) = HalfWingTest(x[1:n], x[n+1:2n], x[2n+1:3n-1], x[3n:4n-2], x[4n-1:end])
winger(x, n) = mean_aerodynamic_chord(wingy(x,n))

cs  = [1.0, 0.6, 0.3]
ts  = [2.0, 2.0, 0.0]
sps = [5.0, 1.0]
dis = [11.3, 60.]
sws = [2.29, 30.]

# winger(x, n) = mean_aerodynamic_chord(HalfWing(chords = x[1:n], twists = x[n+1:2n], spans = x[2n+1:3n-1], dihedrals = x[3n:4n-2], LE_sweeps = x[4n-1:end]))

xs = [cs; ts; sps; dis; sws]
wing = winger(xs, length(cs))

## ForwardDiff (PASSES)
ForwardDiff.gradient(x -> winger(x, length(cs)), xs)

## ReverseDiff (PASSES)
ReverseDiff.gradient(x -> winger(x, length(cs)), xs)

## Zygote (PASSES)
Zygote.gradient(mean_aerodynamic_chord, HalfWingTest(cs, ts, sps, dis, sws))

## FoilerWing
#===================================================================#

# @proto 
struct FoilerWing{T <: Real, N <: AbstractFoil}
    foils  :: Vector{N}
    chords :: Vector{T}
end

FoilerWing(fs :: AbstractVector{N}, cs :: AbstractVector{T}) where {T <: Real, N <: AbstractFoil} = FoilerWing{T,N}(fs, cs)

AeroMDAO.arc_length(fw :: FoilerWing) = sum(arc_length, fw.foils)

# Test
foiler_wing(x1, x2) = arc_length(FoilerWing(fill(Foil(x1), length(x2)), x2))

foilwing = foiler_wing(x_coords, cs)

diff_foiler_wing(x) = foiler_wing(reshape(x[1:end-length(cs)], size(x_coords)), x[end-length(cs):end])
diff_foiler_wing([ x_coords[:]; cs ])

## ForwardDiff (PASSES)
ForwardDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs ])
ForwardDiff.gradient(x -> arc_length(FoilerWing(fill(Foil(x), length(cs)), cs)), x_coords)
ForwardDiff.gradient(x -> arc_length(FoilerWing(fill(Foil(x_coords), length(x)), x)), cs)

## ReverseDiff (PASSES)
ReverseDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs ])
ReverseDiff.gradient(x -> arc_length(FoilerWing(fill(Foil(x), length(cs)), cs)), x_coords)
ReverseDiff.gradient(x -> arc_length(FoilerWing(fill(Foil(x_coords), length(x)), x)), cs)


## Zygote (FAILS)
# ∂f∂x(x, y) = Zygote.gradient(x -> foiler_wing(x, y), x)
# ∂f∂y(x, y) = Zygote.gradient(y -> foiler_wing(x, y), y)
# ∂f∂x(x_coords, cs)
# ∂f∂y(x_coords, cs)

## FoilerWinger
#===================================================================#

# @proto 
struct FoilerWinger{T <: Real, M <: Real, N <: AbstractFoil}
    foils  :: Vector{N}
    chords :: Vector{T}
    twists :: Vector{M}
end

FoilerWing(fs :: AbstractVector{N}, cs :: AbstractVector{T}, ts :: AbstractVector{M}) where {T <: Real, M <: Real, N <: AbstractFoil} = FoilerWing{T,N}(fs, cs, ts)

AeroMDAO.arc_length(fw :: FoilerWinger) = sum(arc_length, fw.foils) + sum(fw.twists)

# Test
foiler_wing(x1, x2, x3) = arc_length(FoilerWinger(fill(Foil(x1), length(x2)), x2, x3))

foilwing = foiler_wing(x_coords, cs, ts)

diff_foiler_wing(x) = foiler_wing(reshape(x[1:end-length(cs)-length(ts)], size(x_coords)), x[end-length(cs)-length(ts):end-length(ts)], x[end-length(ts):end])
diff_foiler_wing([ x_coords[:]; cs; ts ])

## ForwardDiff (PASSES)
ForwardDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs; ts ])
ForwardDiff.gradient(x -> arc_length(FoilerWinger(fill(Foil(x), length(cs)), cs, ts)), x_coords)
ForwardDiff.gradient(x -> arc_length(FoilerWinger(fill(Foil(x_coords), length(x)), x, ts)), cs)
ForwardDiff.gradient(x -> arc_length(FoilerWinger(fill(Foil(x_coords), length(x)), cs, x)), ts)

## ReverseDiff (FAILS)
# ReverseDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs; ts ])
# ReverseDiff.gradient(x -> arc_length(FoilerWinger(fill(Foil(x), length(cs)), cs, ts)), x_coords)
# ReverseDiff.gradient(x -> arc_length(FoilerWinger(fill(Foil(x_coords), length(x)), x, ts)), cs)

## Wing
#===================================================================#

AeroMDAO.arc_length(fw :: HalfWing) = sum(arc_length, fw.foils)
AeroMDAO.arc_length(wing :: Wing) = arc_length(wing.left) + arc_length(wing.right)

winglord(x, n) = HalfWing(foils = fill(Foil(x_coords), n), chords = x[1:n], twists = x[n+1:2n], spans = x[2n+1:3n-1], dihedrals = x[3n:4n-2], LE_sweeps = x[4n-1:end])

wing = winglord(xs, length(cs))

## ForwardDiff (PASSES)
ForwardDiff.gradient(x -> arc_length(winglord(x, length(cs))), xs)
ForwardDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x), length(cs)), cs, ts, sps, dis, sws)), x_coords)
ForwardDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)), x, ts, sps, dis, sws)), cs)
ForwardDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)), cs, x, sps, dis, sws)), ts)
ForwardDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)-1), cs, ts, x, dis, sws)), sps)

## ReverseDiff (PASSES)
ReverseDiff.gradient(x -> arc_length(winglord(x, length(cs))), xs)
ReverseDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x), length(cs)), cs, ts, sps, dis, sws)), x_coords)
ReverseDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)), x, ts, sps, dis, sws)), cs)
ReverseDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)), cs, x, sps, dis, sws)), ts)
# ReverseDiff.gradient(x -> arc_length(HalfWing(fill(Foil(x_coords), length(x)-1), cs, ts, x, dis, sws)), sps)

## VLM
#===================================================================#

using ComponentArrays

function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
    # Evaluate case
    data =  solve_case(aircraft, fs;
                       rho_ref   = ρ,
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print     = print
                      );

    nf_coeffs = sum(nearfield_coefficients(data))
    ff_coeffs = sum(farfield_coefficients(data))

    [ ff_coeffs[1:3]; nf_coeffs[4:6] ]
end

function vlmer(x)
    wing        = winglord(x, length(cs))
    ρ           = 1.225
    ref         = zeros(3)
    V, α, β     = 1.0, 0.0, 0.0
    Ω           = zeros(3)
    fs          = Freestream(V, α, β, Ω)
    S, b, c     = projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
    wing_panels = panel_wing(wing, [20], 10);
    aircraft    = ComponentArray(wing = Horseshoe.(wing_panels...))
    results     = vlm_analysis(aircraft, fs, ρ, ref, S, b, c)
end

vlmer(xs)

## ForwardDiff (PASSES)
ForwardDiff.jacobian(vlmer, xs)