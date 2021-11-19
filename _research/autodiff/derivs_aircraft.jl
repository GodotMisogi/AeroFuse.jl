using Revise
using AeroMDAO
using BenchmarkTools
using StaticArrays
using LinearAlgebra
using ForwardDiff, ReverseDiff, Zygote
using ProtoStructs

## Foil tests
#===================================================================#

# @proto 
struct TestFoil{T <: Real}
    coords :: Matrix{T}
end

# a :: TestFoil + b :: TestFoil = TestFoil(a.coords .+ b.coords)
TestFoil(coords :: AbstractMatrix{T}) where T <: Real = TestFoil{T}(SVector.(coords[:,1], coords[:,2]))

AeroMDAO.arc_length(foil :: TestFoil) = let c = foil.coords; norm(c[2:end] .- c[1:end-1]) end

x_coords = naca4(0,0,1,2)

## ForwardDiff (PASSES)
ForwardDiff.gradient(arc_length ∘ TestFoil, x_coords)

## ReverseDiff (FAILS)
ReverseDiff.gradient(arc_length ∘ TestFoil, x_coords)

## Zygote (FAILS)
Zygote.gradient(arc_length ∘ TestFoil, x_coords) # PASSES FOR MATRIX

## Vector of Foil testss
#===================================================================#

# @proto 
struct TestFoils{T <: Real}
    foils :: Vector{TestFoil{T}}
end

AeroMDAO.arc_length(fs :: TestFoils) = sum(arc_length, fs.foils)

foiler(x) = arc_length(TestFoils(fill(TestFoil(x), 5)))
# airfoils = TestFoils(fill(TestFoil(x_coords), 5))

## ForwardDiff (PASSES)
ForwardDiff.gradient(foiler, x_coords)

## Zygote
# Zygote.gradient(arc_length ∘ TestFoils ∘ (x -> fill(x, 5)), TestFoil(x_coords)) # PASSES
# Zygote.gradient(foiler, x_coords) # FAILS

# x = Zygote.gradient(sum ∘ (x -> reduce(vcat, fill(x, 5))), [1,2,3])

# Great success!

## Wing tests
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

## Zygote (FAILS)
Zygote.gradient(x -> winger(x, length(cs)), xs)

# Great success!

## Composed wing-foil
# @proto 
struct FoilerWing{T <: Real}
    foils  :: Vector{TestFoil{T}}
    chords :: Vector{T}
end

# FoilerWing(fs :: AbstractVector{Foil{T}}, cs :: AbstractVector{T}) where T <: Real = FoilerWing{T}(fs, cs)

# FoilerWing(fs :: AbstractArray{Foil{<: Real}}, cs :: AbstractArray{<: Real}) = FoilerWing(fs, cs)

AeroMDAO.arc_length(fw :: FoilerWing) = sum(arc_length, fw.foils)

# Test
foiler_wing(x1, x2) = arc_length(FoilerWing(fill(TestFoil(x1), length(x2)), x2))

foilwing = foiler_wing(x_coords, cs)

diff_foiler_wing(x) = foiler_wing(reshape(x[1:end-length(cs)], size(x_coords)), x[end-length(cs):end])
diff_foiler_wing([ x_coords[:]; cs ])

## ForwardDiff (PASSES)
ForwardDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs])

## Zygote (FAILS)
∂f∂x(x, y) = Zygote.gradient(x -> foiler_wing(x, y), x)
∂f∂y(x, y) = Zygote.gradient(y -> foiler_wing(x, y), y)
∂f∂x(x_coords, cs)
# ∂f∂y(x_coords, cs)

# Great success!

## Wing
AeroMDAO.arc_length(fw :: HalfWing) = sum(arc_length, fw.foils)
AeroMDAO.arc_length(wing :: Wing) = arc_length(wing.left) + arc_length(wing.right)

winglord(x, n) = Wing(foils = fill(Foil(naca4(2,4,1,2)), n), chords = x[1:n], twists = x[n+1:2n], spans = x[2n+1:3n-1], dihedrals = x[3n:4n-2], LE_sweeps = x[4n-1:end])

winglord(xs, length(cs))

## ForwardDiff (FAILS)
df = ForwardDiff.gradient(arc_length ∘ (x -> winglord(x, length(cs))), xs)

## VLM
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

    nf_coeffs, ff_coeffs = data["Aircraft"][1:2]

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
    aircraft    = Dict("Wing" => wing_panels)
    results     = vlm_analysis(aircraft, fs, ρ, ref, S, b, c)
end

vlmer(xs)

## ForwardDiff
ForwardDiff.jacobian(vlmer, xs)