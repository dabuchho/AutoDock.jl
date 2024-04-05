using Plots

# Complex Numbers & ways to plot them
z = [1 + im, -1 + im]
z = im .* z
scatter(
    z,
    xlims=(-2,2),
    ylims=(-2,2),
    minorgrid=true,
    framestyle= :origin, 
    lab = false,
    size=(500,500)
)

base = Complex[0]
tip = Complex[1 + 1im]

tip = im * tip
quiver(
    real(base), 
    imag(base), 
    quiver=(real(tip), imag(tip)),
    xlims=(-2,2),
    ylims=(-2,2),
    framestyle= :origin,
    size=(500,500)
)