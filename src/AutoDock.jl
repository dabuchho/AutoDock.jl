module AutoDock

# TODO: Only include Package Methods
using Quaternions
using Rotations
using LinearAlgebra
#using CoordinateTransformations # only meant for Vector types. Could implement for BALL Vector3
using StaticArrays
using Intervals: Interval
using ArgCheck
using Distances: sqeuclidean

include("common.jl")

include("random.jl")

include("quaternion.jl")

include("conf.jl")

end
