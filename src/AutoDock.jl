module AutoDock

# TODO: Only include Package Methods
using Quaternions
using Rotations
using LinearAlgebra
using CoordinateTransformations # only meant for Vector types. Could implement for BALL Vector3
using Rotations
using StaticArrays
using Intervals: Interval

include("common.jl")

include("quaternion.jl")

include("conf.jl")

end
