using Zernike
using Zernike: ϵ_fit

# uniform sampling over the unit disk
ρ = range(0.0, 1.0, ϵ_fit)
OPD = stack(2sinc.(5ρ) for i = 1:ϵ_fit; dims = 1)
(; fig, coeffs, latex) = zernike(0, 4)
wait(display(fig))
(; v, fig) = wavefront(OPD, 8)
wait(display(fig))
(; fig) = wavefront(OPD, [(0, 2), (0, 4), (0, 6), (0, 8)])
wait(display(fig))
(; fig) = transform(v, 0.2)
wait(display(fig))
(; fig) = transform(v, 0.2, 0.17 + 0.17im)
wait(display(fig))
(; fig) = transform(v, 0.2, 0.0im, 0.0, (0.83, 4.18))
wait(display(fig))
(; fig) = transform(v, 0.2, 0.17 + 0.17im, float(π), (0.83, 4.18))
wait(display(fig))
