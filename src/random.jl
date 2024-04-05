# This does NOT work with Integer bound boxes!
function random_inside_box(corner1::Vector{T}, corner2::Vector{T}) where T
    tmp::Vector{T} = zeros(T, 3)
    for i in eachindex(tmp)
        tmp[i] = rand(T) * (corner2[i] - corner1[i]) + corner1[i]
    end # corner1[i] < corner2[i]
    tmp
end

 # sphere of radius 1 centered in the origin
function random_inside_sphere(type::T = Float64) where T
    tmp = rand(type, 3) .* 2 .- 1
    if sqr(tmp) < 1     # sqrt() not required, since radius is 1
        return tmp
    else
        random_inside_sphere(type)
    end
end


# consider adding a function rand(lower_bound, upper_bound)