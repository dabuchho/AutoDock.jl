@testitem "Quaternions" begin
    using Quaternions
    using LinearAlgebra
    using Rotations
    # type tests    
    for T in [Float32, Float64]
        q = random_quaternion(T)    # normalized
        @test typeof(q) == Quaternion{T}

        @test typeof(quaternion_norm_sqr(q)) == T
        @test typeof(quaternion_normalize(q)) == Quaternion{T}

        @test quaternion_is_normalized(q) == true
        q_not_normalized = generate_quaternion(randn(T),randn(T),randn(T),randn(T))
        @test quaternion_is_normalized(q_not_normalized) == false

        v = [randn(T), randn(T), randn(T)] * 10 # quaternion_from_angle normalizes v
        ϑ = randn(T) * 10 # quaternion_from_angle normalizes ϑ
        @test typeof(quaternion_from_angle(v, ϑ)) == Quaternion{T}
        @test typeof(quaternion_from_angle(v))    == Quaternion{T}
        @test typeof(vector_from_quaternion(q))   == Vector{T}
        @test typeof(rotation_from_quaternion(q)) == QuatRotation{T}
        @test typeof(random_quaternion())         == QuaternionF64
        @test typeof(random_quaternion(T))        == Quaternion{T}
        @test typeof(quaternion_increment(q, v))  == Quaternion{T}

        q2 = random_quaternion(T)
        @test typeof(quaternion_difference(q,q2)) == Vector{T}

        # correctness
        IQ = Quaternion{T}(1,0,0,0)  # identity Quaternion
        v₀ = zeros(T, 3)             # zero Vector
        ϑ  = rand(T)                 # random angle
        
        # TODO: merge type tests with correctness tests
        @test qt_identity(T) == IQ
        #= Consider writing a test that compares "random" Quaternions (seeded)
           generate_quaternion() == random_quaternion()
        =#
        r = rand(T, 4)
        q = generate_quaternion(r[1], r[2], r[3], r[4])

        @test quaternion_norm_sqr(q) == sum(r.*r)
        q_normalized = quaternion_normalize(q)
        @test length(q_normalized) == 1
        quaternion_is_normalized(q_normalized) == true  # testing this twice
        
        @test quaternion_from_angle(v₀)    == IQ
        @test quaternion_from_angle(v₀, ϑ) == IQ
        
        # TODO: predict outcome (for a given set of parameters for example)
        v = rand(T, 3)
        @test quaternion_from_angle(v, ϑ) != IQ
        @test quaternion_from_angle(v)    != IQ
        
        v_small = T[eps(T) * T(0.1), 0, 0]    # norm(v_small) < ε
        @test quaternion_from_angle(v_small) == IQ

        @test vector_from_quaternion(IQ)           == v₀  
        @test vector_from_quaternion(q_normalized) != v₀ 

        @test rotation_from_quaternion(IQ) == Matrix(I, 3, 3)
        @test rotation_from_quaternion(q)  != Matrix(I, 3, 3)
        # random_quaternion() 
        # quaternion_increment() 
        # quaternion_dfference()
    end
end

using Quaternions
using LinearAlgebra
using Rotations

# I = Quaternion{Float32}(1,0,0,0)
# q = generate_quaternion(1,0,0,0)
# QuatRotation{Float32}(q) == Matrix(I, 3, 3)