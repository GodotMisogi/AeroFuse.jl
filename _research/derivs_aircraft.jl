using Revise
using Rotations
using AeroMDAO
using BenchmarkTools
using ForwardDiff
using StaticArrays

## Wing tests
struct HalfWingTest{T <: Real}
    chords    :: Vector{T}
    twists    :: Vector{T}
    spans     :: Vector{T}
    dihedrals :: Vector{T}
    sweeps    :: Vector{T} 
    HalfWingTest(chords :: AbstractVector{T}, twists :: AbstractVector{T}, spans :: AbstractVector{T}, dihedrals :: AbstractVector{T}, sweeps :: AbstractVector{T}) where T <: Real = new{T}(chords, -deg2rad.(twists), spans, deg2rad.(dihedrals), deg2rad.(sweeps))
end

cs  = [1.0, 0.6]
ts  = [2.0, 2.0]
sps = [5.0]
dis = [11.3]
sws = [2.29]

wing = HalfWingTest(cs, ts, sps, dis, sws)
winger(x, n) = HalfWingTest(x[1:n], x[n+1:2n], x[2n+1:3n-1], x[3n:4n-2], x[])
ForwardDiff.gradient(mean_aerodynamic_center ∘ winger, xs)

## Old tests
TrapezoidalWing(b, δ, Λ, λ, c_root, τ_root, τ_tip) = WingSection(span       = b,
																 dihedral 	= δ,
																 sweep_LE   = Λ,
																 taper      = λ,
																 root_chord = c_root,
																 root_twist = τ_root,
																 tip_twist  = τ_tip)

xs = [4.0, 0.0, 15.0, 0.4, 2.0, 0.0, -2.0]
wing  	 = TrapezoidalWing(xs...)
wing_mac = mean_aerodynamic_center(wing)
print_info(wing, "Wing")

## Differentiation tests
winger(x) = TrapezoidalWing(x...)
ForwardDiff.gradient(mean_aerodynamic_center ∘ winger, xs)

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
	wing 		= winger(x)
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