using Revise
using AeroMDAO
using BenchmarkTools
using ForwardDiff, ReverseDiff
using StaticArrays
using LinearAlgebra

## Foil tests 
struct FoilTest{T <: Real}
	foils :: Vector{Foil{T}}
end

arc_length(foil :: Foil) = let c = foil.coords; norm(c[2:end] .- c[1:end-1]) end

arc_length(fs :: FoilTest) = sum(arc_length, fs.foils)

foiler(x) = arc_length(FoilTest(fill(Foil(x, "NACA 0012"), 5)))

x_coords = let coords = naca4(0,0,1,2); [ getindex.(coords, 1) getindex.(coords, 2) ] end
airfoils = FoilTest(fill(Foil(x_coords, "NACA 0012"), 5))

foiler(x_coords)
ReverseDiff.gradient(foiler, x_coords)

# Great success!

## Wing tests
struct HalfWingTest{T <: Real}
    chords    :: Vector{T}
    twists    :: Vector{T}
    spans     :: Vector{T}
    dihedrals :: Vector{T}
    sweeps    :: Vector{T} 
    HalfWingTest(chords :: AbstractVector{T}, twists :: AbstractVector{T}, spans :: AbstractVector{T}, dihedrals :: AbstractVector{T}, sweeps :: AbstractVector{T}) where T <: Real = new{T}(chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps))
end

AeroMDAO.AircraftGeometry.mean_aerodynamic_chord(root_chord, taper_ratio) = (2/3) * root_chord * (1 + taper_ratio + taper_ratio^2)/(1 + taper_ratio)
section_macs(wing :: HalfWingTest) = @views mean_aerodynamic_chord.(wing.chords[1:end-1], fwddiv(wing.chords))
section_projected_areas(wing :: HalfWingTest) = wing.spans .* fwdsum(wing.chords) / 2

"""
    mean_aerodynamic_chord(half_wing :: HalfWing)
    
Compute the mean aerodynamic chord of a `HalfWing`.
"""
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

# winger(x, n) = mean_aerodynamic_chord(HalfWing(chords = x[1:n], twists = x[n+1:2n], spans = x[2n+1:3n-1], dihedrals = x[3n:4n-2], sweep_LEs = x[4n-1:end]))

xs = [cs; ts; sps; dis; sws]
wing = winger(xs, length(cs)) 

ForwardDiff.gradient(x -> winger(x, length(cs)), xs)

# Great success!

## Composed wing-foil
struct FoilerWing{T <: Real}
    foils  :: Vector{Foil{T}}
    chords :: Vector{T}
    FoilerWing(fs :: AbstractVector{Foil{<: Real}}, cs :: AbstractVector{T}) where T <: Real = new{T}(fs, cs)
end

FoilerWing(fs :: AbstractArray{Foil{<: Real}}, cs :: AbstractArray{<: Real}) = FoilerWing(fs, cs)

arc_length(fw :: FoilerWing) = sum(arc_length, fw.foils)

# Test
foiler_wing(x1, x2) = arc_length(FoilerWing(fill(Foil(x1), length(x2)), x2))

foilwing = foiler_wing(x_coords, cs)

diff_foiler_wing(x) = foiler_wing(reshape(x[1:end-length(cs)], size(x_coords)), x[end-length(cs):end])
diff_foiler_wing([ x_coords[:]; cs ])
∂f∂x(x, y) = ForwardDiff.gradient(x -> foiler_wing(x, y), x)
∂f∂y(x, y) = ForwardDiff.gradient(y -> foiler_wing(x, y), y)
ForwardDiff.gradient(diff_foiler_wing, [ x_coords[:]; cs])
∂f∂x(x_coords, cs)
∂f∂y(x_coords, cs)

# Great success!

## Wing
winglord(x, n) = Wing(chords = x[1:n], twists = x[n+1:2n], spans = x[2n+1:3n-1], dihedrals = x[3n:4n-2], sweep_LEs = x[4n-1:end])


## VLM
function vlm_analysis(aircraft, fs, ρ, ref, S, b, c, print = false)
	# Evaluate case
	data = 	solve_case(aircraft, fs; 
                       rho_ref   = ρ, 
                       r_ref     = ref,
                       area_ref  = S,
                       span_ref  = b,
                       chord_ref = c,
                       print 	 = print
                      );
		
	nf_coeffs, ff_coeffs, CFs, CMs, horseshoe_panels, normals, horseshoes, Γs = data["Aircraft"]
			   
	[ ff_coeffs[1:3]; nf_coeffs[4:6] ]
end

function vlmer(x)
	wing 		= winger(x, length(cs))
	wing_pos    = [0., 0., 0.]
	ρ 			= 1.225
	ref     	= zeros(3)
	V, α, β 	= 1.0, 0.0, 0.0
	Ω 			= zeros(3)
	fs 	    	= Freestream(V, α, β, Ω)
	S, b, c 	= projected_area(wing), span(wing), mean_aerodynamic_chord(wing)
	wing_panels = panel_wing(wing, [20], 10,
                          	 position = wing_pos);
	aircraft 	= Dict("Wing" => wing_panels)
	results 	= vlm_analysis(aircraft, fs, ρ, ref, S, b, c)
end

vlmer(xs)

## ForwardDiff
ForwardDiff.jacobian(vlmer, [5.0, 0.0, 10.0, 0.4, 2.0, 0.0, -2.0])