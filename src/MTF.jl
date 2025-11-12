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

const mtf_map = Dict(
    :x => (1:div(d_max, 2), 1),
    :y => (1, 1:div(d_max, 2))
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
