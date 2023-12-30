using Zernike
using Zernike: Z, W, P, ϵ_fit

# uniform sampling over the unit disk
ρ = range(0.0, 1.0, ϵ_fit)
OPD = stack(2sinc.(5ρ) for i = 1:ϵ_fit; dims = 1)
(; fig, coeffs, latex) = Z(0, 4)
wait(display(fig))
recap_a, v, metrics1, fig, axis, plot = W(OPD, 8)
wait(display(fig))
(; fig) = W(OPD, [(0, 2), (0, 4), (0, 6), (0, 8)])
wait(display(fig))
recap_b, v2, metrics2, fig = P(v, 0.2)
wait(display(fig))
(; fig) = P(v, 0.2, 0.17 + 0.17im)
wait(display(fig))
(; fig) = P(v, 0.2, 0.0im, 0.0, (0.83, 4.18))
wait(display(fig))
(; fig) = P(v, 0.2, 0.17 + 0.17im, float(π), (0.83, 4.18))
wait(display(fig))
