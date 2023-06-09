using Zernike

# uniform sampling over the unit disk
ρ = range(0.0, 1.0, 21)
OPD = stack(2sinc.(5ρ) for i = 1:21; dims = 1)
fig, coeffs, latex = Z(0, 4)
wait(fig.scene.current_screens[])
a, v, metrics1, fig = W(OPD, 8)
wait(fig.scene.current_screens[])
(; fig) = W(OPD, [(0, 2), (0, 4), (0, 6), (0, 8)])
wait(fig.scene.current_screens[])
b, v2, metrics2, fig = P(v, 0.2)
wait(fig.scene.current_screens[])
(; fig) = P(v, 0.2, 0.17 + 0.17im)
wait(fig.scene.current_screens[])
(; fig) = P(v, 0.2, 0.0im, 0.0, (0.83, 4.18))
wait(fig.scene.current_screens[])
(; fig) = P(v, 0.2, 0.17 + 0.17im, float(π), (0.83, 4.18))
wait(fig.scene.current_screens[])
