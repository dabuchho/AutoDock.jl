using Rotations
using Quaternions
using BenchmarkTools


function rot_by_matrix(v::Vector, α::Real)
    RotX(α) * v
end


function rot_by_quatrot(v::Vector, α::Real)
    QuatRotation(cos(α/2),sin(α/2),0,0) * v
end


function quat_from_axisangle(axis::AbstractVector, theta::Real)
    if length(axis) != 3
        error("Must be a 3-vector")
    end
    s, c = sincos(theta / 2)
    axis = normalize(axis)
    return Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
end



α, β, γ = 1.2, -0.8, 0.1;
x = [1,1,1]

xm = rot_by_matrix(x, π)
xq = rot_by_quatrot(x, π)


@btime rot_by_matrix(x, α)

@btime rot_by_quatrot(x, α)


using Plots
gr()

scatter([x[1]], [x[2]], [x[3]], 
    xlab = "x", 
    ylab = "y", 
    zlab= "z",
    xlims=(-2,2),
    ylims=(-2,2),
    zlims=(-2,2)
)
scatter!([xm.x], [xm.y], [xm.z])
scatter!([xq.x], [xq.y], [xq.z])


qr = QuatRotation(1,2,3,4)
inv = [0,0,0]
out = [1,1,1]
qr * inv + out
