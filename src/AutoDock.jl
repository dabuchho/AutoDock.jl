module AutoDock

# TODO: Only include Package Methods
using Quaternions
using Rotations
using LinearAlgebra
#using CoordinateTransformations # only meant for Vector types. Could implement for BALL Vector3
using StaticArrays
using Intervals: Interval
using ArgCheck

include("common.jl")

include("quaternion.jl")
#export quaternion_norm_sqr
include("conf.jl")

end
