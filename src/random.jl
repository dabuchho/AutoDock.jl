# This does NOT work with Integer bound boxes!
function random_inside_box(corner1::Vector{T}, corner2::Vector{T}) where T
    tmp::Vector{T} = zeros(T, 3)
    for i in eachindex(tmp)
        tmp[i] = rand(T) * (corner2[i] - corner1[i]) + corner1[i]
    end # corner1[i] < corner2[i]
    tmp
end

# consider adding a function rand(lower_bound, upper_bound)