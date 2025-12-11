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

const mtf_map = Dict(
    :x => (1:div(d_max, 2), 1),
    :y => (1, 1:div(d_max, 2))
)

const psf_map = Dict(
    :x => (:, div(d_max, 2) + 1),
    :y => (div(d_max, 2) + 1, :)
)

function mtf_plot(MTF::FloatMat)
    (; colormap, fontsize) = plotconfig
    ξx = ξy = range(-1.0, 1.0, size(MTF, 1))
    figure = (; plotconfig.size)
    axis = (
        type = Axis3,
        title = "Modulation Transfer Function",
        zlabel = L"MTF",
        xlabel = L"\xi_x",
        ylabel = L"\xi_y",
        titlesize = fontsize,
        xlabelsize = fontsize,
        ylabelsize = fontsize,
        zlabelsize = fontsize,
    )
    with_theme(plotconfig.theme) do
        surface(ξx, ξy, MTF; figure, axis, colormap)
    end
end

function mtf_plot(MTF::FloatMat, x_or_y::Symbol)
    MTF = fftshift(MTF)
    ξ = range(0.0, 1.0, div(size(MTF, 1), 2))
    s = "Modulation Transfer Function"
    title = x_or_y === :x ? "Sagittal $s" : "Tangential $s"
    xlabel = "Normalized frequency"
    ylabel = "MTF"
    plot = lines(ξ, view(MTF, mtf_map[x_or_y]...); axis = (; title, xlabel, ylabel))
    DataInspector()
    plot
end

function mtf_plot!(MTF::FloatMat, x_or_y::Symbol)
    MTF = fftshift(MTF)
    ξ = range(0.0, 1.0, div(size(MTF, 1), 2))
    lines!(ξ, view(MTF, mtf_map[x_or_y]...))
end

function psf_plot(PSF::FloatMat)
    (; colormap, fontsize) = plotconfig
    x = y = range(-1.0, 1.0, size(PSF, 1))
    figure = (; plotconfig.size)
    axis = (
        type = Axis3,
        title = "Point Spread Function",
        zlabel = L"PSF",
        xlabel = L"x",
        ylabel = L"y",
        titlesize = fontsize,
        xlabelsize = fontsize,
        ylabelsize = fontsize,
        zlabelsize = fontsize,
    )
    with_theme(plotconfig.theme) do
        surface(x, y, PSF; figure, axis, colormap)
    end
end

function psf_plot(PSF::FloatMat, x_or_y::Symbol)
    c = range(-1.0, 1.0, size(PSF, 1))
    s = "Point Spread Function"
    title = x_or_y === :x ? "Sagittal $s" : "Tangential $s"
    xlabel = string(x_or_y)
    ylabel = "PSF"
    plot = lines(c, view(PSF, psf_map[x_or_y]...); axis = (; title, xlabel, ylabel))
    DataInspector()
    plot
end

function psf_plot!(PSF::FloatMat, x_or_y::Symbol)
    c = range(-1.0, 1.0, size(PSF, 1))
    lines!(c, view(PSF, psf_map[x_or_y]...))
end
