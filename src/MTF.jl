using FFTW

function PSF(ΔW::Wavefront; s = 50.0)
    a = 2.0 * s
    x = y = range(-a, a, d_max)
    P = [hypot(x, y) ≤ 1.0 ? ei2pi(-ΔW(complex(x, y))) : 0.0im for x ∈ x, y ∈ y]
    psf = P |> fft |> fftshift .|> abs2
    psf .*= metrics(ΔW).strehl / maximum(psf)
    return psf
end

OTF(ΔW::Wavefront) = PSF(ΔW; s = 1.0) |> fft |> fftshift

function MTF(ΔW::Wavefront)
    mtf = abs.(OTF(ΔW))
    mtf ./= maximum(mtf)
    return mtf
end
