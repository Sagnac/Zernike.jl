import GLMakie: GLFW.GetPrimaryMonitor, GLFW.GetMonitorContentScale,
                MonitorProperties, activate!

mutable struct PlotConfig
    size::Tuple{Float64, Float64}
    fontsize::Float64
    colormap::Symbol
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
        focus_on_show = true,
    )
end

const plotconfig = PlotConfig(default_config()...)

function refresh_monitor()
    plotconfig.size, plotconfig.fontsize = get_monitor_properties()
end

function reset_config()
    defaults = default_config()
    for name in fieldnames(PlotConfig)
        setfield!(plotconfig, name, defaults[name])
    end
end

function setproperty!(plotconfig::PlotConfig, name::Symbol, value)
    if name == :reset && value
        reset_config()
    elseif name == :resize && value
        refresh_monitor()
    else
        setfield!(plotconfig, name, convert(fieldtype(PlotConfig, name), value))
    end
end

propertynames(plotconfig::PlotConfig) = fieldnames(PlotConfig)..., :reset, :resize

function ZPlot(ρ::FloatVec, θ::FloatVec, Zp::FloatMat;
               high_order = false, window = "ZernikePlot", plot = window)
    (; size, fontsize, colormap, focus_on_show) = plotconfig
    axis3attributes = (
        title = plot,
        titlesize = fontsize,
        xlabel = L"\rho_x",
        ylabel = L"\rho_y",
        zlabel = L"Z",
        xlabelsize = fontsize,
        ylabelsize = fontsize,
        zlabelsize = fontsize,
        azimuth = 0.275π,
        protrusions = 80,
    )
    ρᵪ = [ρⱼ * cos(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
    ρᵧ = [ρⱼ * sin(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
    fig = Figure(; size)
    axis3 = Axis3(fig[1,1]; axis3attributes...)
    if high_order
        @. Zp = sign(Zp) * log10(abs(Zp * log(10)) + 1)
        axis3.title[] = "Log transform of $plot"
    end
    surface!(axis3, ρᵪ, ρᵧ, Zp; shading = NoShading, colormap)
    # hacky way to produce a top-down heatmap-style view without generating
    # another plot with a different set of data
    # accomplished by adding a toggle which changes the perspective on demand
    zproperties = (
        :zlabelvisible,
        :zgridvisible,
        :zticksvisible,
        :zticklabelsvisible,
        :zspinesvisible,
    )
    label = Label(fig, "Pupil view"; fontsize = 0.76fontsize)
    pupil = Toggle(fig; active = true)
    fig[1,2] = grid!([label pupil]; tellheight = false, valign = :bottom)
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
    register_interaction!(axis3, :resize) do event::MouseEvent, axis::Axis3
        if event.type === MouseEventTypes.rightclick &&
           ispressed(axis.scene, Keyboard.left_control)
            resize_to_layout!(fig)
            return Consume(true)
        else
            return Consume(false)
        end
    end
    activate!(; title = window, focus_on_show)
    display(fig)
    return fig
end
