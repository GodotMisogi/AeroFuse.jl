## Packages
using AeroFuse
using Plots
using Accessors
using LinearAlgebra
using StaticArrays

## Testing
panel = Panel3D(
    [1.0, -1.0, 0.0], 
    [0.0, -1.0, 0.0], 
    [0.0, 0.0,  0.0], 
    [1.0, 0.0,  0.0]
)

##
quadrilateral_doublet_potential(1.0, panel, [0,0,0])

##
quadrilateral_source_potential(1.0, panel, [0,0,-1e-15])

## Wing analysis case
wing = WingSection(
    root_foil  = naca4(0,0,1,2),
    tip_foil   = naca4(0,0,1,2),
    taper      = 1.0,
    aspect     = 6.0,
    dihedral   = 0.0,
    sweep      = 0.0,
    symmetry   = true
)

wing_mesh = WingMesh(wing, [12], 6, span_spacing = Uniform());

##
@time surf_mesh = surface_coordinates(wing_mesh, wing_mesh.num_span, 10)
@time surf_pans = make_panels(surf_mesh);
@time wake_pans = map(surf_pans[end,:]) do p
    temp_p = setproperties(p,
        p1 = p.p2,
        p4 = p.p3
    )
    fin1_p = @set temp_p.p2[1] = 10wing.chords[1] 
    fin2_p = @set fin1_p.p3[1] = 10wing.chords[1]
end

## 
using Plots

plot(surf_pans, 
    camera = (30,30),
    aspect_ratio = 1, 
    zlim = (-0.5, 0.5) .* span(wing), 
    lc = :cornflowerblue,
)

plot!([wake_pans;; ], lc = :grey)


## Build aircraft
ac = ComponentVector(
    wing = surf_pans,
    wake = wake_pans,
)

Nc = size(surf_pans, 1)
Ns = size(surf_pans, 2)
Nw = Ns

## Aerodynamic Influence Coefficient matrix
AIC_ff = [ if pan_i === pan_j; 0.5 else quadrilateral_doublet_potential(1., pan_i, collocation_point(pan_j)) end for pan_i in ac.wing[:], pan_j in ac.wing[:] ]

##
AIC_fw = [ quadrilateral_doublet_potential(1., pan_i, collocation_point(pan_j)) for pan_i in ac.wing[:], pan_j in ac.wake[:] ]


##
AIC_ff[:,1] -= AIC_fw[:,1]
for i in 1:Nw
    AIC_ff[:,i+Nc] -= AIC_fw[:,i]
end

##
AIC_k = zeros(Nw, Nc * Ns)

for i = 1:Nw
    for j = 1:(Nw * Nc)
        if i == j
            AIC_k[i,j] = 1
        elseif j == i + Nc
            AIC_k[i,j] = -1
        end
    end
end

AIC_k

##
AIC = @views [
    AIC_ff AIC_fw;
    AIC_k   -I(Nw)
]

##
BIC_ff = [ quadrilateral_source_potential(1., pan_i, collocation_point(pan_j)) for pan_i in ac.wing[:], pan_j in ac.wing[:] ]

## Boundary condition
V = SVector(1,0,0)
RHS = -[
    [ dot(V, collocation_point(p_j)) for p_j in ac.wing[:] ];
    zeros(Nw)
]

## Solve
phi = AIC \ RHS

reshape(phi[1:(Nc * Ns)], Nc, Ns)

##
