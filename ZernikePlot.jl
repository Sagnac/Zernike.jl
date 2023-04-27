const ρᵪ = [ρⱼ * cos(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]
const ρᵧ = [ρⱼ * sin(θᵢ) for θᵢ ∈ θ, ρⱼ ∈ ρ]

function ZPlot(Zp; titles...)

    resolution = (1400, 1000)
    textsize = 35

    axis3attributes = (
        title = L"%$(titles[:plot])",
        titlesize = textsize,
        xlabel = L"\rho_x",
        ylabel = L"\rho_y",
        zlabel = L"Z",
        xlabelsize = textsize,
        ylabelsize = textsize,
        zlabelsize = textsize,
        azimuth = 0.275π,
        protrusions = 80,
    )

    fig = Figure(; resolution)
    axis3 = Axis3(fig[1,1]; axis3attributes...)
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

    label = Label(fig, "Pupil view"; textsize = 0.76textsize)
    pupil = Toggle(fig)
    fig[1,2] = grid!([label pupil]; tellheight = false, valign = :bottom)

    on(pupil.active) do active
        for property in zproperties
            getfield(axis3, property)[] = !active
        end
        axis3.azimuth = active ? -π/2 : 0.275π
        axis3.elevation = active ? π/2 : π/8
        axis3.ylabeloffset = active ? 90 : 40
        axis3.xlabeloffset = active ? 50 : 40
    end

    set_window_config!(title = titles[:window], focus_on_show = true)
    display(fig)

    return

end
