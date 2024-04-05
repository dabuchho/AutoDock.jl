@testitem "Quaternions" begin
    using Quaternions
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
        # ...
    end
end