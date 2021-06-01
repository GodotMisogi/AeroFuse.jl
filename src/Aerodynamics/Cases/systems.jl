struct VLMAnalysis{T <: Real}
	aircraft 	:: Dict{String, NTuple{2, Matrix{Panel3D{T}}}} 
	freestream 	:: Freestream{T}
	r_ref 		:: SVector{3,T}
	rho_ref 	:: T
	area_ref 	:: T
	chord_ref 	:: T
	span_ref 	:: T
	name 		:: String
end

VLMAnalysis(components :: Dict{String, NTuple{2, Matrix{Panel3D{T}}}}, freestream :: Freestream; rho_ref = 1.225, r_ref = zeros(3), area_ref = 1, chord_ref = 1, span_ref = 1, name = "Aircraft") where T <: Real = VLMAnalysis{T}(components, freestream, r_ref, rho_ref, area_ref, chord_ref, span_ref, name)

mutable struct VLMSystem{T <: Real}
    AIC         :: Matrix{T}
    RHS         :: Vector{T}
    nf_coeffs   :: SVector{9,T}
    ff_coeffs   :: SVector{3,T}
    CFs         :: Matrix{SVector{3,T}}
    CMs         :: Matrix{SVector{3,T}}
    horseshoes  :: Matrix{Horseshoe{T}}
    normals     :: Matrix{SVector{3,T}}
    Γs          :: Matrix{T}
end

function solve_system(vlm_anal :: VLMAnalysis) 
    horseshoe_panels, camber_panels = values(vlm_anal.aircraft)
    normals = panel_normal.(camber_panels)

    solve_case(horseshoe_panels, normals, vlm_anal.freestream, vlm_anal.rho_ref, vlm_anal.r_ref)
end