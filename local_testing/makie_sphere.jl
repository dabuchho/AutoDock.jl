###############################################################################
#    This script is a simple test in preparation for Quaternion operations    #
###############################################################################
using GLMakie, Random, Colors, LinearAlgebra
using GeometryBasics
using Makie

# Sphere with arrow
with_theme(theme_black()) do
    fig = Figure(size = (600,600))
    ax = LScene(fig[1,1], show_axis=true) # [1,1] centers plot
    sphere = Sphere(Point3f(0), 1) # Not sure yet what these parameters do
    # Sphere
    meshscatter!(
        ax,
        Point3f(0,0,0),  # Center
        sphere = sphere, 
        markersize = 0.025, # Radius
        transparency = true)
    # Vector
    arrows!(
        ax,
        [Point3f(0,0,0)],  # Base
        [Point3f(0,0,0.05)],  # Head
        arrowsize = Vec3f(0.02, 0.02, 0.02),
        color = :red,
        arrowcolor = :red)
    fig
end
GLMakie.closeall() # close any open screen