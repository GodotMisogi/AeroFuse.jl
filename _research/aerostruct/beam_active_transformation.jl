# Assemble RHS
function assemble_fem_dynamics(Fx, Fy, Fz, Mx, My, Mz)
    px = [ Fx; Mx ]                                 # Fx Mx concatenated
    py = (collect ∘ Iterators.flatten ∘ zip)(Fy, My) # Fy My interspersed
    pz = (collect ∘ Iterators.flatten ∘ zip)(Fz, Mz) # Fz Mz interspersed

    F  = [ py; pz; px ]
end

function assemble_fem_dynamics(loads)
    px = [ loads[1,:]; loads[4,:] ] # Fx Mx concatenated
    py = loads[[2,5],:][:]          # Fy My interspersed
    pz = loads[[3,6],:][:]          # Fz Mz interspersed
    F  = [ py; pz; px ]
end


## "FEM" setup
stiffness_matrix = tube_stiffness_matrix(aluminum, tubes)

## Loads
load_vector = assemble_fem_dynamics(Fx, Fy, Fz, Mx, My, Mz)


## Solve system(s)
xs = stiffness_matrix \ load_vector

function compute_displacements(Δs, n)
    # Assembly
    n1 = 2n                 # dim(Fy + My)
    n2 = n1 + 2n            # dim(Fz + Mz)
    n3 = n2 +  n            # dim(Fx)
    n4 = n3 +  n            # dim(Mx)

    dy = @view Δs[1:2:n1]
    θy = @view Δs[2:2:n1]
    dz = @view Δs[n1+1:2:n2]
    θz = @view Δs[n1+2:2:n2]
    dx = @view Δs[n2+1:n3]
    θx = @view Δs[n3+1:n4]

    dx, θx, dy, θy, dz, θz
end

dx, θx, dy, θy, dz, θz = compute_displacements(xs, length(fem_pts))


## Axis system transformation matrices
#==========================================================================================#

# Change of basis for representations in the first frame
orthogonal_change_of_basis(second, third)        = second * third * second'
orthogonal_change_of_basis(first, second, third) = orthogonal_change_of_basis(second, third') * first

# Compute local beam axes
ss = @. normalize(fem_pts[2:end] - fem_pts[1:end-1])
ns = fill(SVector(1,0,0), length(fem_pts) - 1)
# normalize((plan_coords[end,2:end] - plan_coords[1,1:end-1]) × (plan_coords[1,2:end] - plan_coords[end,1:end-1]))
cs = @. ss × ns

# Section axes
# bound_coords = coordinates(wing)
# bound_pts    = (1 - fem_w) * bound_coords[1,:] + fem_w * bound_coords[end,:]
# ss = @. normalize(bound_pts[2:end] - bound_pts[1:end-1])
# ns = @. normalize((bound_coords[end,2:end] - bound_coords[1,1:end-1]) × (bound_coords[1,2:end] - bound_coords[end,1:end-1]))
# cs = @. ss × ns

## Active transformations of forces and moments (vs. passive transformation of stiffness matrix)
global_axis     = I(3)            # Global axis system for panels in VLM
beam_local_axis = [ 0. -1.  0. ;
                   -1.  0.  0. ;
                    0.  0. -1. ]  # Orthogonal local axis system for beam in FEM

# WTF array of local coordinate systems
wtf = zeros(3, 3, size(ss)...) # ((x,y,z), (c,s,n), chordwise, spanwise)
for inds in CartesianIndices(wtf)
    l,k,i  = inds.I
    wtf[l,1,i] = cs[i][l]
    wtf[l,2,i] = ss[i][l]
    wtf[l,3,i] = ns[i][l]
end
wtf

## Array comprehension
mats = @views [ orthogonal_change_of_basis(global_axis, wtf[:,:,inds], beam_local_axis) for inds in CartesianIndices(ss) ]

# Einstein summation
# @einsum mats_eins[l,k,i,j] := local_axes[m,l,i,j] * global_axis[m,k]
## Transform forces
fem_forces  = [ mats .* pt_forces[1:end-1] ; [mats[end] * pt_forces[end] ] ]
fem_moments = [ mats .* pt_moments[1:end-1]; [mats[end] * pt_moments[end]] ]
# fem_loads   = [ reduce(hcat, fem_forces); reduce(hcat, fem_moments) ]

# fem_forces  = [ Ref(mats[1]) .* pt_forces[1:span_num-1] ; zero_vec; Ref(mats[end]) .* pt_forces[span_num+1:end]  ]
# fem_moments = [ Ref(mats[1]) .* pt_moments[1:span_num-1]; zero_vec; Ref(mats[end]) .* pt_moments[span_num+1:end] ]

## Concatenate forces and moments into loads array
fem_loads   = [ reduce(hcat, fem_forces); reduce(hcat, fem_moments) ]