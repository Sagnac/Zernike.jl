import GLMakie: GLFW.GetPrimaryMonitor, GLFW.GetMonitorContentScale,
                MonitorProperties, activate!, Screen

mutable struct PlotConfig
    size::Tuple{Float64, Float64}
    fontsize::Float64
    colormap::Symbol
    theme::Makie.Attributes
    focus_on_show::Bool
end

function get_monitor_properties()
    monitor = GetPrimaryMonitor()
    monitor_properties = MonitorProperties(monitor)
    monitor_scale = GetMonitorContentScale(monitor)
    dpi_scale = sum(monitor_scale) / length(monitor_scale)
    (; height) = monitor_properties.videomode
    size = (0.72height, 0.6height) ./ dpi_scale
    fontsize = 0.13 * monitor_properties.dpi[1] / dpi_scale
    return (; size, fontsize)
end

function default_config()
    return (;
        get_monitor_properties()...,
        colormap = :oslo,
        theme = theme_black(),
        focus_on_show = true,
    )
end

const plotconfig = PlotConfig(default_config()...)

function resize!(plotconfig::PlotConfig)
    plotconfig.size, plotconfig.fontsize = get_monitor_properties()
end

function reset!(plotconfig::PlotConfig; only_theme = false)
    set_theme!(; GLMakie = (; title = "Makie", focus_on_show = false))
    only_theme && return plotconfig
    defaults = default_config()
    for name in fieldnames(PlotConfig)
        setfield!(plotconfig, name, defaults[name])
    end
    return plotconfig
end

const window_title = "ZernikePlot"

const shading = NoShading

const zproperties = (
    :zlabelvisible,
    :zgridvisible,
    :zticksvisible,
    :zticklabelsvisible,
    :zspinesvisible,
)

function zernikeplot!(axis, φ; m = 10, n = 10, finesse = finesse,
                      high_order = false, colormap = plotconfig.colormap, kwargs...)
    φ = convert(Observable, φ)
    ρ, θ = polar(m, n; finesse)
    x, y = polar_mat(ρ, θ)
    z = @lift($φ.(ρ', θ))
    if high_order
        zv = z[]
        ln10 = log(10.0)
        @. z[] = sign(zv) * log10(abs(zv * ln10) + 1.0)
    end
    surface!(axis, x, y, z; shading, colormap, kwargs...)
end

function zernikeplot!(axis, ρ, θ, φ::FloatMat; kwargs...)
    ρ = vec(ρ)
    θ = vec(θ)
    x, y = polar_mat(ρ, θ)
    colormap = get(kwargs, :colormap, plotconfig.colormap)
    surface!(axis, x, y, φ; shading, colormap, kwargs...)
end

function zernikeplot!(axis, ρ, θ, φ; kwargs...)
    ρ = vec(ρ)
    θ = vec(θ)
    zernikeplot!(axis, ρ, θ, φ.(ρ', θ); kwargs...)
end

function _zplot(args...; plot_title = window_title,
                size = plotconfig.size,
                fontsize = plotconfig.fontsize,
                kwargs...)
    if get(kwargs, :high_order, false)
        plot_title = LaTeXString("Log transform of " * plot_title)
    end
    axis3attributes = (
        title = plot_title,
        titlesize = fontsize,
        xlabel = L"x",
        ylabel = L"y",
        zlabel = L"Z",
        xlabelsize = fontsize,
        ylabelsize = fontsize,
        zlabelsize = fontsize,
        azimuth = 0.275π,
        protrusions = 80,
    )
    fig = Figure(; size)
    axis3::Axis3 = Axis3(fig[1,1]; axis3attributes...)
    plot::SurfacePlot = zernikeplot!(axis3, args...; kwargs...)
    on(_ -> reset_limits!(axis3), plot[3])
    # hacky way to produce a top-down heatmap-style view without generating
    # another plot with a different set of data
    # accomplished by adding a toggle which changes the perspective on demand
    pgrid = GridLayout(fig[1,2]; tellheight = false, valign = :bottom)
    Label(pgrid[1,1], "Pupil view"; fontsize = 0.76fontsize)
    pupil::Toggle = Toggle(pgrid[2,1]; active = true)
    on(pupil.active; update = true) do active
        for property in zproperties
            getproperty(axis3, property)[] = !active
        end
        axis3.azimuth = active ? -π/2 : 0.275π
        axis3.elevation = active ? π/2 : π/8
        axis3.ylabeloffset = active ? 90 : 40
        axis3.xlabeloffset = active ? 50 : 40
    end
    colsize!(fig.layout, 1, Aspect(1, 1.0))
    resize_to_layout!(fig)
    onmouserightup(_ -> resize_to_layout!(fig), addmouseevents!(fig.scene))
    return FigureAxisPlot(fig, axis3, plot)
end

function zplot(args...;
    window_title = window_title,
    focus_on_show = plotconfig.focus_on_show,
    theme = plotconfig.theme,
    kwargs...
)
    screen_config = (; title = window_title, focus_on_show)
    if Makie.current_backend() != GLMakie
        GLMakie.activate!(; screen_config...)
    else
        set_theme!(; GLMakie = screen_config)
    end
    with_theme(() -> _zplot(args...; kwargs...), theme)
end

(φ::Phase)() = zplot(φ)

(φ::Phase)(::Type{T}) where T <: Screen = display(T(), φ())

(φ::Observable{<:Phase})() = φ[]()
