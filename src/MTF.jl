using FFTW

function PSF(ΔW::Wavefront)
    x = y = range(-2.0, 2.0, d_max)
    P = [hypot(x, y) ≤ 1.0 ? ei2pi(-ΔW(complex(x, y))) : 0.0im for x ∈ x, y ∈ y]
    P |> fft |> fftshift .|> abs2
end

OTF(ΔW::Wavefront) = PSF(ΔW) |> fft |> fftshift

function MTF(ΔW::Wavefront)
    mtf = abs.(OTF(ΔW))
    mtf / maximum(mtf)
end
