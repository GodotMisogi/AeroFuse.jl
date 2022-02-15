abstract type AbstractFreestream end

"""
    Freestream(α, β, Ω)
    Freestream(α, β, Ω_x, Ω_y, Ω_z)
    Freestream(U, Ω)
    Freestream(; 
                 alpha = 0., 
                 beta  = 0., 
                 omega = [0,0,0]
               )

Define freestream conditions with angle of attack ``α`` (degrees), sideslip angle ``β`` (degrees), and (quasi-steady) rotation vector ``Ω`` for a vortex lattice analysis.

Alternatively, provide the velocity vector ``U``, which is normalized to determine the angles.
"""
struct Freestream{T} <: AbstractFreestream
    alpha :: T
    beta  :: T
    omega :: SVector{3,T}
end

Base.broadcastable(fs :: Freestream) = Ref(fs)

Freestream(α_deg, β_deg, Ω) = let T = promote_type(eltype(α_deg), eltype(β_deg), eltype(Ω)); Freestream{T}(deg2rad(α_deg), deg2rad(β_deg), Ω) end

Freestream(α_deg, β_deg, Ω_x, Ω_y, Ω_z) = Freestream(deg2rad(α_deg), deg2rad(β_deg), SVector(Ω_x, Ω_y, Ω_z))

Freestream(U, Ω) = let (V, α, β) = cartesian_to_freestream(normalize(U)); Freestream{T}(V, α, β, Ω) end

Freestream(; alpha = 0., beta = 0., omega = zeros(3)) = Freestream(alpha, beta, omega)

function Base.show(io :: IO, fs :: Freestream)
    println(io, "Freestream: ")
    for fname in fieldnames(typeof(fs))
        println(io, "    ", fname, " = ", getfield(fs, fname))
    end
end

"""
    velocity(freestream :: Freestream)

Compute the velocity of a `Freestream`.
"""
velocity(fs :: Freestream) = freestream_to_cartesian(1., fs.alpha, fs.beta)

"""
    body_frame_velocity(freestream :: Freestream)

Compute the velocity of Freestream in the body reference frame.
"""
body_frame_velocity(fs :: Freestream) = -velocity(fs)
