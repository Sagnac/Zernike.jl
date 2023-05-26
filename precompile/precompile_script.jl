using Zernike

# uniform sampling over the unit disk
ρ = repeat(range(0, 1, 21); outer = 21)

θ = repeat(range(0, 2π, 21); inner = 21)

OPD = 2sinc.(5ρ)

fig, coeffs, latex = Z(0, 4)

wait(fig.scene.current_screens[])

a, v, metrics, fig = W(ρ, θ, OPD, 8)

wait(fig.scene.current_screens[])

b, v2, metrics, fig = S(0.2, v)

wait(fig.scene.current_screens[])
