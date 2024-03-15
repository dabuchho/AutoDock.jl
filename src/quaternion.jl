# Further Reading
# https://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
# https://arxiv.org/pdf/math/0701759.pdf


function quaternion_norm_sqr(q::Quaternion)
    q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3    
end


# same as normalize(q::Quaternion) but faster
function quaternion_normalize(q::Quaternion)
    s = quaternion_norm_sqr(q)
    a = sqrt(s)
    q *= 1/a # slightly faster than q /= a
end


# triple redundant
function quaternion_is_normalized(q::Quaternion)
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
function quaternion_from_angle(axis::AbstractVector, ϑ::Number)
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    ϑ = normalize_angle(ϑ)
    s, c = sincos(ϑ / 2)
    axis = normalize(axis)
    Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
end

function quaternion_from_angle(v::AbstractVector)
    ϑ = norm(v, 2)              # 2-Norm of v
    if ϑ > eps(Float64)
        axis = (1 / ϑ) * v
        return quaternion_from_angle(axis, ϑ)
    else
        Quaternion(1, 0, 0, 0)  # Identity Quaternion
    end
end

# named quaternion_to_angle in AutoDock
function vector_from_quaternion(q::Quaternion)
    if quaternion_is_normalized(q)
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
    else
        error("Quaternion has to be normalized")
    end
end


# named quaternion_to_r3 in AutoDock
function rotation_from_quaternion(q::Quaternion)
    QuatRotation(q)
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


# named random_orientation in AutoDock
function random_quaternion()
    q = Quaternion(     # normal distributed (μ=0, σ=1)
        randn(),
        randn(),
        randn(),
        randn()
    )
    if norm(q) > eps(Float64)   
        q = quaternion_normalize(q)
    else                # ≊ (0, 0, 0, 0), so generate a different Quaternion
        random_quaternion()
    end 
end


function quaternion_increment(q::Quaternion, v::AbstractVector)
    if !quaternion_is_normalized(q)
        error("Quaternion has to be normalized")
    end
    q = quaternion_from_angle(v) * q
    quaternion_normalize(q)
end


# rotation that needs to be applied to convert a to b
function quaternion_difference(b::Quaternion, a::Quaternion)
    if !quaternion_is_normalized(a) || !quaternion_is_normalized(b)
        error("Quaternions have to be normalized")
    end
    tmp = b
    tmp /= a
    vector_from_quaternion(tmp)
end