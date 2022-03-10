## INVISCID SYSTEM SETUPS
#============================================#

## Doublet-source solution and edge velocities
#============================================#

function solve_inviscid_doublets(panels, wakes, u)
    # Inviscid solution
    all_panels      = [ panels; wakes ]
    foil_doublets   = doublet_matrix(panels, panels)
    all_sources     = source_matrix(panels, all_panels)
    A_inv           = foil_doublets^(-1)
    P               = -A_inv * all_sources

    # Mass defect block computation
    all_lengths = panel_length.(all_panels)         #  N + N_w
    B           = source_matrix(wakes, all_panels)  # (N_w, N + N_w)
    A_w         = doublet_matrix(wakes, panels)     # (N_w, N)
    dKds        = -midpair_map(-, defect_block(P, A_w, B)) ./ all_lengths
    D           = dKds ./ all_lengths'

    # Inviscid velocities
    φs          = A_inv * boundary_vector(panels, u, (0., 0.))
    q0s, Δrs    = tangential_velocities(panels, φs, u, false)
    U_invs      = [ q0s; Ref(u) .+ tangent_vector.(wakes) ]

    U_invs, D, all_lengths
end

## Linear-vortex solution and edge velocities
#============================================#

function solve_inviscid_vortices(panels, wakes, u)
    # Inviscid solution
    all_panels      = [ panels; wakes ]
    foil_vortices   = vortex_influence_matrix(panels)
    all_sources     = source_influence_matrix(panels, all_panels)
    A_inv           = foil_vortices^(-1)
    P               = -A_inv * all_sources

    # Mass defect block computation
    all_lengths     = panel_length.(all_panels)                 # Foil and Wake: N + N_w
    B               = linear_source_neumann_matrix(wakes, all_panels)
    A_w             = linear_vortex_neumann_matrix(wakes, panels)
    C               = defect_block(P, A_w, B)[:,1:end-1]
    D               = C ./ all_lengths'

    # Inviscid velocities
    U_inv_f     = A_inv * neumann_boundary_condition(panels, u) # Foil: N + 1
    U_inv_w     = neumann_boundary_condition(wakes, u)          # Wake: N_w + 1
    U_invs      = [ U_inv_f[1:end-1]; U_inv_w ]                 # Foil and wake: N + N_w + 2

    U_invs, D, all_lengths
end

## DEFECT FORMULATION
#============================================#

defect_block(P, Aw, B)         = [ P; Aw * P + B ]
edge_velocities(ms, D, U_invs) = U_invs + D * ms
