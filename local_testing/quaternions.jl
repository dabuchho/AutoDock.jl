using Quaternions
using Plots, Plots.Measures

# Basics
k = quat(0, 0, 0, 1)

j = quat(0, 0, 1, 0)

i = j*k

i^2 == j^2 == k^2 == i*j*k == -1

1 + i + k + j

test = quat(1, 2, 3, 4)
print(test.s, test.v1, test.v2, test.v3) # How to access each postition

# Operations on Quaternions
function rotate_vector(q::Quaternion, u::AbstractVector)
    if length(u) != 3
        error("Must be a 3-vector")
    end
    q_u = Quaternion(0, u[1], u[2], u[3])
    q_v = q*q_u*conj(q)
    return [imag_part(q_v)...]
end

struct Q{T}
    _q::Quaternion{T}
    function Q{T}(q::T,w::T,e::T,r::T) where T <:Number
        _q = Quaternion(q,w,e,r)
        new(_q)
    end
end

Q(q::T,w::T,e::T,r::T) where T<:Number = Q{T}(q,w,e,r)

function Q(q::T,w::T,e::T,r::T) where T
    Q(q,w,e,r)
end

q = Q(1,1,1,1)
q._q

################################ PLOTS #########################################
# Simple 3d Quiver Plot
x = [0]
y = [0]
z = [0]

u = [1]
v = [1]
w = [1]
quiver(x,y,z,quiver=(u,v,w))

# Sphere Quiver Plot
ϕs = range(-π, π, length = 50)
θs = range(0, π, length = 25)
θqs = range(1, π - 1, length = 25)
x = vec([sin(θ) * cos(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
y = vec([sin(θ) * sin(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
z = vec([cos(θ) for (ϕ, θ) = Iterators.product(ϕs, θs)])
u = 0.1 * vec([sin(θ) * cos(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
v = 0.1 * vec([sin(θ) * sin(ϕ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
w = 0.1 * vec([cos(θ) for (ϕ, θ) = Iterators.product(ϕs, θqs)])
quiver(x, y, z, quiver = (u, v, w))


using Quaternions

function quat(q::T,w::T,e::T,r::T) where T
    Quaternion(q,w,e,r)
end

a::Float32 = 1.0
b::Float32 = 2.0
c::Float32 = 3.0
d::Float32 = 4.0
q = quat(a,b,c,d)



function quaternion_norm_sqr(q::Quaternion{T}) where T
    q.s*q.s + q.v1*q.v1 + q.v2*q.v2 + q.v3*q.v3    
end

typeof(quaternion_norm_sqr(q))

q2 = QuaternionF64(1.0,1.0,1.0,2.0)
typeof(quaternion_norm_sqr(q2))

