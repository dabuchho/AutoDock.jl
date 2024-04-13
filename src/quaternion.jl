# Further Reading
# https://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
# https://arxiv.org/pdf/math/0701759.pdf

export
    qt_identity,
    generate_quaternion,
    quaternion_norm_sqr, 
    quaternion_normalize, 
    quaternion_is_normalized,
    quaternion_from_angle, 
    vector_from_quaternion, 
    rotation_from_quaternion,
    random_quaternion, 
    quaternion_increment, 
    quaternion_difference

    
function qt_identity(type::T) where T
    Quaternion{type}(1)
end

# only used for testing
function generate_quaternion(a::T,b::T,c::T,d::T)::Quaternion{T} where T
    Quaternion{T}(a,b,c,d)
end


function quaternion_norm_sqr(q::Quaternion{T})::T where T
    q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3    
end


# same as normalize(q::Quaternion) but faster
function quaternion_normalize(q::Quaternion{T})::Quaternion{T} where T
    s::T = quaternion_norm_sqr(q)
    a::T = sqrt(s)
    q *= 1/a    # slightly faster than q /= a
end


# triple redundant
function quaternion_is_normalized(q::Quaternion{T})::Bool where T
    quaternion_norm_sqr(q) ≈ 1 && norm(q) ≈ 1 && abs(q) ≈ 1 
end


# Operations on Quaternions
# from Documentation Quaternions.jl
# function rotate_vector(q::Quaternion, u::AbstractVector)
#     if length(u) != 3
#         error("Must be a 3-vector")
#     end
#     q_u = Quaternion(0, u[1], u[2], u[3])
#     q_v = q*q_u*conj(q)
#     return [imag_part(q_v)...]
# end


# from Documentation Quaternions.jl (named angle_to_quaternion in AutoDock)
# ϑ must be a Number, because π is typed as Irrational
function quaternion_from_angle(axis::Vector{T}, ϑ::T)::Quaternion{T} where T
    @assert length(axis) == 3
    ϑ = normalize_angle(ϑ)
    s, c = sincos(ϑ / 2)
    axis = normalize(axis)
    Quaternion{T}(c, s*axis[1], s*axis[2], s*axis[3])
end

a = [2.0f0,3.0f0,1.0f0]
sig = 4.7f0


function quaternion_from_angle(v::Vector{T})::Quaternion{T} where T
    ϑ = norm(v, 2)              # 2-Norm of v
    if ϑ > eps(T)
        axis = (1 / ϑ) * v
        return quaternion_from_angle(axis, ϑ)
    else
        Quaternion{T}(1.0, 0.0, 0.0, 0.0)  # Identity Quaternion
    end
end

# named quaternion_to_angle in AutoDock
function vector_from_quaternion(q::Quaternion{T})::Vector{T} where T
    @assert quaternion_is_normalized(q)
    c = q.s
    if c > -1 && c < 1
        ϑ = 2*acos(c)  # cos⁻¹(c) ∈ [0, +π]
        if ϑ > π
            ϑ -= 2π    # ϑ ∈ [-π, +π]
        end
        v = [q.v1, q.v2, q.v3]
        s = sin(ϑ/2)
        if abs(s) < eps(Float64)
            return zeros(3)
        end
        return v *= ϑ / s
    else
        return zeros(3)
    end
end


# named quaternion_to_r3 in AutoDock
function rotation_from_quaternion(q::Quaternion{T})::QuatRotation{T} where T
    QuatRotation{T}(q)
end


#= 
named quaternion_to_r3 in AutoDock
TODO: QuatRotation(::Quaternion) already implemented this. But automatically
normalizes Quaternions, if they are not normalized before.
=#        
# function rotation_from_quaternion(q::Quaternion)
#     if !quaternion_is_normalized(q)
#         error("Quaternion has to be normalized")
#     end
#     # TODO: shorten this
#     a = q.s
#     b = q.v1
#     c = q.v2
#     d = q.v3

#     aa = a*a
#     ab = a*b
#     ac = a*c
#     ad = a*d
#     bb = b*b
#     bc = b*c
#     bd = b*d
#     cc = c*c
#     cd = c*d
#     dd = d*d

#     if !isapprox(aa+bb+cc+dd, 1)
#         error("Quaternion has to be of 2-Norm")
#     end

#     # https://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
#     # 4- Quaternions and rotations
    
#     # QuatRotations are Row-major ordered (and immutable?)
#     rot = QuatRotation(    
#         aa+bb-cc-dd,    # rot[1][1]
#         2*(ad+bc),      # rot[2][1]
#         2*(-ac+bd),     # ...
        
#         2*(-ad+bc),     # rot[1][2]
#         aa-bb+cc-dd,    # ...
#         2*(ab+cd),
        
#         2*(ac+bd),
#         2*(-ab+cd),
#         aa-bb-cc+dd
#     )
# end


# named random_orientation in AutoDock # TODO: add seed
function random_quaternion(type = Float64)::Quaternion{type}
    q = Quaternion{type}(     # normal distributed (μ=0, σ=1)
        randn(type),
        randn(type),
        randn(type),
        randn(type)
    )
    if norm(q) > eps(type)   
        q = quaternion_normalize(q)
    else                # ≊ (0, 0, 0, 0), so generate a different Quaternion
        random_quaternion()
    end 
end                     # this is not very pretty


function quaternion_increment(q::Quaternion{T}, v::Vector{T})::Quaternion{T} where T
    @assert quaternion_is_normalized(q)
    q = quaternion_from_angle(v) * q
    quaternion_normalize(q)
end


# rotation that needs to be applied to convert a to b
function quaternion_difference(
    b::Quaternion{T}, 
    a::Quaternion{T}
)::Vector{T} where T

    @assert quaternion_is_normalized(b)
    @assert quaternion_is_normalized(a)
    tmp = b
    tmp /= a
    vector_from_quaternion(tmp)
end