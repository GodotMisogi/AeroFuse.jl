abstract type AbstractFreestream end

struct Freestream{T} <: AbstractFreestream
    alpha :: T
    beta  :: T
    omega :: SVector{3,T}
end


Freestream(α_deg, β_deg, Ω) = let T = promote_type(eltype(α_deg), eltype(β_deg)); Freestream{T}(deg2rad(α_deg), deg2rad(β_deg), Ω) end

Freestream(α_deg, β_deg, Ω_x, Ω_y, Ω_z) = Freestream(deg2rad(α_deg), deg2rad(β_deg), SVector(Ω_x, Ω_y, Ω_z))

"""
    Freestream(U, Ω)

A Freestream flow in Cartesian coordinates with vector ``U`` and quasi-steady rotation vector ``\\Omega``.
"""
Freestream(U, Ω) = let (V, α, β) = cartesian_to_freestream(U); Freestream{T}(V, α, β, Ω) end

"""
    Freestream(V, α, β, Ω)
    
A Freestream flow in spherical polar coordinates with magnitude ``V``, angle-of-attack ``α``, side-slip angle ``β``, and a quasi-steady rotation vector ``\\Omega``.
"""
Freestream(; alpha = 0., beta = 0., omega = zeros(3)) = Freestream(alpha, beta, omega)


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