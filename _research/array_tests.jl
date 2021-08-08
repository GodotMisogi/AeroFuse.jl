
# Matrix versions
ac_fem   = reshape(permutedims(reduce(hcat, ac_mine)), (size(ac_mine)..., 3))
fem_arr  = (1 - fem_w) * mesh[1,:,:] + fem_w * mesh[end,:,:]

# Moments
M_in  = [ let xs = (ac_fem[i,:,:] - fem_arr[1:end-1,:]); 
              @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
              for i in eachindex(ac_fem[:,1,1]) ]
M_out = [ let xs = (ac_fem[i,:,:] - fem_arr[2:end,:]); 
              @. SVector(xs[:,1], xs[:,2], xs[:,3]) × forces[i,:] / 2 end
              for i in eachindex(ac_fem[:,1,1]) ]

M_ins  = sum(M_in)
M_outs = sum(M_out)

# Load setups
fem_loads  = zeros(6, length(pt_forces))
fem_forces = reduce(hcat, half_forces)
fem_loads[1:3,1:end-1] = fem_forces
fem_loads[1:3,2:end]  += fem_forces

fem_M_ins  = reduce(hcat, M_ins)
fem_M_outs = reduce(hcat, M_outs)
fem_loads[4:end,1:end-1] = fem_M_ins
fem_loads[4:end,2:end]  += fem_M_outs

# Constrain the plane of symmetry and print
fem_loads[:,span_num+1] .= 0.
fem_loads
 
# Transform loads to principal axes
fem_force_trans  = [ dircos[1] * fem_loads[1:3,1:span_num]   dircos[2] * fem_loads[1:3,span_num+1:end]   ]
fem_moment_trans = [ dircos[1] * fem_loads[4:end,1:span_num] dircos[2] * fem_loads[4:end,span_num+1:end] ]
fem_loads_trans  = [ fem_force_trans  ; 
                     fem_moment_trans ]
fem_loads_trans

## Splitting for Wing case (useless now?)
#==========================================================================================#

function middle_index(x :: AbstractArray)
    n = length(x)
    if n % 2 == 0
        Int(n / 2)
    else
        ceil(Int, n / 2)
    end
end

middle(x :: AbstractVector) = x[middle_index(x)]

zero_vec = [SVector(0,0,0.)]

## Testing permutations and transformations of stiffness matrices
#==========================================================================================#

using Einsum
using SparseArrays
using UnicodePlots

# Beam properties
E     = 85e9  # Elastic modulus, N/m²
G     = 25e9  # Shear modulus, N/m²
σ_max = 350e6 # Yield stress with factor of safety 2.5, N/m²
ρ     = 1.6e3 # Density, kg/m³
ν     = 0.3   # Poisson ratio (UNUSED FOR NOW)
R     = 1e-2  # Outer radius, m
t     = 8e-3  # Thickness, m
# Ls    =  [ 0.9610428525060863
#            0.8669692821235053
#            0.6880307178764946
#            0.4417429131718409
#            0.1522142343220727 ]

# Create material and tubes
aluminum = Material(E, G, σ_max, ρ)
tubes    = Tube.(Ref(aluminum), Ls, R, t)

# Testing blocks
Ks = [ tube_stiffness_matrix(aluminum, [tube]) for tube in tubes ]
Ks[1]

## Testing simultaneous permutations of rows and columns:
# 1. (Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2): [9,1,5,11,6,2,10,3,7,12,8,4]
# 2. (Fx1, Fx2, Mx1, Mx2, Fy1, My1, Fy2, My2, Fz1, Mz1, Fz2, Mz2): [9,10,11,12,1,2,3,4,5,6,7,8]
inds    = [9,1,5,11,6,2,10,3,7,12,8,4]
perm_Ks = [ K[inds,:][:,inds] for K in Ks ] 

# Axis transformation of stiffness matrices
axis_trans = [ 0  1  0; 
               0  0 -1; 
              -1  0  0]

# Tensoring for force and moments following 1.
K_trans = kron(I(4), axis_trans)

K_tran = repeat(K_trans, 1, 1, length(Ls))

perm_K = zeros(12, 12, length(Ls))
for i in 1:length(Ls)
    perm_K[:,:,i] .= perm_Ks[i]
end

## Using Einstein summation convention for multiplication: D = TᵗKT

@einsum D[j,k,i] := K_tran[l,j,i] * perm_K[l,m,i] * K_tran[m,k,i]


# Constraints
cons = [span_num]

stiffness = build_stiffness_matrix(D, cons)

UnicodePlots.spy(stiffness)

## Concatenate forces and moments into loads array
fem_loads   = [ reduce(hcat, pt_forces); reduce(hcat, pt_moments) ]

x = solve_cantilever_beam(D, fem_loads, cons)

## Array tests
vlm_arr = reshape(reduce(hcat, bound_leg_center.(horsies)), (3, size(horsies)...))
vlm_fs  = reshape(reduce(hcat, vlm_forces), (3, size(horsies)...))
fem_arr = reduce(hcat, fem_pts)

mommy(vlm_arr, fem_arr, vlm_fs) = 
    [ (vlm_arr[:,j,i] - fem_arr[:,i]) × vlm_fs[:,j,i] for j in eachindex(vlm_fs[1,:,1]), i in eachindex(fem_arr[1,:]) ] 

cuck = sum(mommy(vlm_arr, fem_arr[:,1:end-1], vlm_fs / 2), dims = 1)[:]


    # dx_wtf      = reshape(reduce(hcat, dx_mesh), (size(dx_mesh)..., 3))
    # shape       = size(moment_arms)
    # mommy_wtf   = reshape(reduce(hcat, moment_arms), (shape[1], shape[2], 3))
    
    # @einsum wtf[k,l,i] := T[l,i,j] * mommy_wtf[k,l,j]
