import GLMakie: GLFW.GetPrimaryMonitor, MonitorProperties

function ZPlot(ρ, θ, Zp; high_order = false, window = "ZernikePlot", plot = window)
    monitor_properties = MonitorProperties(GetPrimaryMonitor())
    (; height) = monitor_properties.videomode
    resolution = (0.85height, height/2)
    fontsize = 0.13 * monitor_properties.dpi[1]
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
    fig = Figure(; resolution)
    axis3 = Axis3(fig[1,1]; axis3attributes...)
    if high_order
        @. Zp = sign(Zp) * log10(abs(Zp * log(10)) + 1)
        axis3.title[] = "Log transform of $plot"
    end
    surface!(axis3, ρᵪ, ρᵧ, Zp; colormap = :oslo)
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
            getfield(axis3, property)[] = !active
        end
        axis3.azimuth = active ? -π/2 : 0.275π
        axis3.elevation = active ? π/2 : π/8
        axis3.ylabeloffset = active ? 90 : 40
        axis3.xlabeloffset = active ? 50 : 40
    end
    set_window_config!(title = window, focus_on_show = true)
    display(fig)
    return fig
end
