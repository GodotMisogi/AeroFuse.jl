# VLM analysis
function vlm_analysis(Γs, new_horsies, other_horsies, ρ, U, Ω, syms)
    # Compute component forces for structural residual
    @timeit "Get Circulations" new_Γs = view.(Ref(Γs), syms) 
    @timeit "Get Aerodynamic Centers" new_acs = map(horsies -> bound_leg_center.(horsies), new_horsies)

    # Allocating steps

    # Collect all horseshoes
    @timeit "All Horseshoes" all_horsies = [ mapreduce(vec, vcat, new_horsies); vec(other_horsies) ]

    # Compute modified horseshoe forces
    @timeit "New Forces" new_forces = map((Γ_comp, hs_comp) -> surface_forces(hs_comp, Γ_comp, all_horsies, Γs, U, Ω, ρ), new_Γs, new_horsies)

    # Compute unmodified horseshoe forces
    @timeit "Get Symbols" other_syms = filter(sym -> !(sym ∈ syms), keys(Γs))
    @timeit "Get Other Circulations" other_Γs = mapreduce(sym -> vec(view(Γs, sym)), vcat, other_syms)
    @timeit "Other Forces" other_forces = surface_forces(vec(other_horsies), other_Γs, all_horsies, Γs, U, Ω, ρ)

    # Get all forces
    @timeit "All Forces" vlm_forces = [ mapreduce(vec, vcat, new_forces); vec(other_forces) ]

    # Non-allocating steps, fails with ForwardDiff for obvious reasons.
    # [ setproperty!(all_horsies, sym, horsies) for (sym, horsies) in zip(syms, new_horsies) ]
    # vlm_forces = surface_forces(all_horsies, Γs, all_horsies, Γs, U, Ω, ρ)
    # new_forces = view.(Ref(vlm_forces), syms)

    return new_forces, vlm_forces, new_acs, all_horsies
end
