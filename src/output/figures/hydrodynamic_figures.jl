"""
    hydrodynamic_figures.jl

Hydrodynamic Eulerian figures for TethysDrift's --viz path.

These are TethysDrift's local versions of the ocean figure functions. The two packages
do not import each other (SPLIT_PLAN.md section 1), so TethysDrift renders hydrodynamic
fields from the store the ocean half wrote, using its own `extract_hydrodynamic_dataset`
reader.
"""

using CairoMakie
import CairoMakie: record

"""
    _finite_colorrange(A; fallback = (0.0, 1.0)) -> Tuple{Float64, Float64}

Compute a finite `(lo, hi)` colorrange from `A`, ignoring `NaN` and `Inf` entries
(such as masked land cells) and expanding degenerate constant fields by ±0.5.
"""
function _finite_colorrange(A; fallback::Tuple{Float64, Float64} = (0.0, 1.0))
    vals = filter(isfinite, vec(A))
    isempty(vals) && return fallback
    lo, hi = Float64.(extrema(vals))
    return lo == hi ? (lo - 0.5, hi + 0.5) : (lo, hi)
end

@inline _surface_depth_index(depths::AbstractVector, depth_level::Union{Nothing, Int}) =
    isnothing(depth_level) ? argmax(depths) : clamp(depth_level, 1, length(depths))

"""
    plot_hydrodynamic_advection(
        hydro::NamedTuple;
        depth_level::Int = 1,
        time_index::Int = 1,
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        title::AbstractString = "Surface Current Velocity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
        quiver_stride::Int = 3,
        colormap::Symbol = :turbo
    ) -> Figure

Plot horizontal current velocity vectors over current speed magnitude at a depth level.
"""
function plot_hydrodynamic_advection(
    hydro::NamedTuple;
    depth_level::Union{Nothing, Int} = nothing,
    time_index::Int = 1,
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    title::AbstractString = "Surface Current Velocity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_advection.png",
    quiver_stride::Int = 3,
    colormap::Symbol = :turbo
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    u_k = Float64.(hydro.u[:, :, k] .* 100.0)  # m/s -> cm/s
    v_k = Float64.(hydro.v[:, :, k] .* 100.0)
    spd_k = Float64.(hydro.speed[:, :, k] .* 100.0)

    fig = Figure(size = (980, 720), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = "$title\n[Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )

    max_spd = let vals = filter(isfinite, spd_k)
        isempty(vals) ? 1.0 : max(1.0, maximum(vals))
    end
    hm = heatmap!(ax, lons, lats, spd_k, colormap = colormap;
                   colorrange = (0, max_spd))
    Colorbar(fig[1, 2], hm, label = "Current Speed (cm/s)")

    # Quiver overlay
    stride = max(1, quiver_stride)
    ix = 1:stride:length(lons)
    iy = 1:stride:length(lats)
    uq = u_k[ix, iy]
    vq = v_k[ix, iy]
    arrows2d!(
        ax,
        lons[ix],
        lats[iy],
        uq,
        vq;
        tipwidth = 8,
        tiplength = 8,
        color = :white
    )

    # Bathymetry contours
    if !isnothing(bathymetry_data) && hasproperty(bathymetry_data, :elevation)
        contour!(
            ax,
            bathymetry_data.lon,
            bathymetry_data.lat,
            bathymetry_data.elevation;
            levels = [-2000, -1000, -500, -200],
            color = :black,
            linewidth = 0.6,
            linestyle = :dash
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_tracers(
        hydro::NamedTuple;
        depth_level::Int = 1,
        time_index::Int = 1,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png",
        colormap_temp::Symbol = :thermal,
        colormap_salt::Symbol = :haline
    ) -> Figure

Plot surface temperature and salinity side by side.
"""
function plot_hydrodynamic_tracers(
    hydro::NamedTuple;
    depth_level::Union{Nothing, Int} = nothing,
    time_index::Int = 1,
    title::AbstractString = "",
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_tracers.png",
    colormap_temp::Symbol = :thermal,
    colormap_salt::Symbol = :haline
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    T_k = Float64.(hydro.temperature[:, :, k])
    S_k = Float64.(hydro.salinity[:, :, k])

    fig = Figure(size = (1200, 480), fontsize = 12)

    ax1 = Axis(fig[1, 1],
        title = "Surface Temperature (°C)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )
    hm1 = heatmap!(ax1, lons, lats, T_k, colormap = colormap_temp,
                   colorrange = _finite_colorrange(T_k; fallback = (0.0, 15.0)))
    Colorbar(fig[1, 2], hm1, label = "°C")

    ax2 = Axis(fig[1, 3],
        title = "Surface Salinity (PSU)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )
    hm2 = heatmap!(ax2, lons, lats, S_k, colormap = colormap_salt,
                   colorrange = _finite_colorrange(S_k; fallback = (30.0, 35.0)))
    Colorbar(fig[1, 4], hm2, label = "PSU")

    fig[0, :] = Label(fig, "Hydrodynamic Tracers - $(isempty(title) ? "" : title)\\n[Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        fontsize = 14, font = :bold)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_stratification(
        hydro::NamedTuple;
        depth_level::Int = 1,
        time_index::Int = 1,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_stratification.png"
    ) -> Figure

Plot buoyancy frequency squared N² at a depth level.
"""
function plot_hydrodynamic_stratification(
    hydro::NamedTuple;
    depth_level::Union{Nothing, Int} = nothing,
    time_index::Int = 1,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_stratification.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    N2_k = Float64.(hydro.stratification[:, :, k] .* 1e4)  # convert to 10^-4 s^-2

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Buoyancy Frequency N² — [Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )
    hm = heatmap!(ax, lons, lats, N2_k, colormap = :ice;
                   colorrange = _finite_colorrange(N2_k; fallback = (0.0, 5.0)))
    Colorbar(fig[1, 2], hm, label = "10⁻⁴ s⁻²")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_diffusion(
        hydro::NamedTuple;
        depth_level::Int = 1,
        time_index::Int = 1,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_diffusion.png"
    ) -> Figure

Plot vertical eddy diffusivity κ_v at a depth level.
"""
function plot_hydrodynamic_diffusion(
    hydro::NamedTuple;
    depth_level::Union{Nothing, Int} = nothing,
    time_index::Int = 1,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_diffusion.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    κ_k = Float64.(hydro.diffusion[:, :, k])

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Vertical Eddy Diffusivity κ_v — [Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )
    κ_log = log10.(clamp.(κ_k, 1e-8, 1e-2))
    hm = heatmap!(ax, lons, lats, κ_log, colormap = :turbid;
                   colorrange = _finite_colorrange(κ_log; fallback = (-8.0, -2.0)))
    Colorbar(fig[1, 2], hm, label = "log10(κ_v) [m²/s]")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_section(
        hydro::NamedTuple;
        variable::Symbol = :temperature,
        coordinate::Real = 44.0,
        section_type::Symbol = :lat,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_section.png"
    ) -> Figure

Plot a vertical cross-section of temperature or salinity along a latitude or longitude line.
"""
function plot_hydrodynamic_section(
    hydro::NamedTuple;
    variable::Symbol = :temperature,
    coordinate::Real = 44.0,
    section_type::Symbol = :lat,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_section.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths

    field = variable === :temperature ? hydro.temperature : hydro.salinity

    if section_type === :lat
        j = argmin(abs.(lats .- coordinate))
        sec = Float64.(field[:, j, :])
        x_axis = lons
        x_label = "Longitude (°E)"
        title = "Temperature Cross-Section at Latitude $(round(coordinate, digits=1))°N"
    else
        i = argmin(abs.(lons .- coordinate))
        sec = Float64.(field[i, :, :])
        x_axis = lats
        x_label = "Latitude (°N)"
        title = "Temperature Cross-Section at Longitude $(round(coordinate, digits=1))°W"
    end

    # Mask below bathymetry
    bathy = hydro.bathymetry
    for ix in 1:length(x_axis), iz in 1:length(depths)
        if section_type === :lat
            b_val = bathy[ix, j]
        else
            b_val = bathy[i, ix]
        end
        if !isnan(b_val) && depths[iz] < b_val
            sec[ix, iz] = NaN
        end
    end

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = title,
        xlabel = x_label,
        ylabel = "Depth (m)"
    )
    hm = heatmap!(ax, x_axis, depths, sec;
                  colormap = :thermal,
                  colorrange = _finite_colorrange(sec; fallback = (0.0, 15.0)))
    Colorbar(fig[1, 2], hm, label = variable === :temperature ? "°C" : "PSU")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_timeseries(
        hydro::NamedTuple;
        station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        variable::Symbol = :temperature,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
    ) -> Figure

Plot Eulerian time series at a representative station.
"""
function plot_hydrodynamic_timeseries(
    hydro::NamedTuple;
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    variable::Symbol = :temperature,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths

    # Find nearest station
    if isnothing(station)
        i_st = length(lons) ÷ 2
        j_st = length(lats) ÷ 2
    else
        i_st = argmin(abs.(lons .- station[1]))
        j_st = argmin(abs.(lats .- station[2]))
    end

    k_surf = 1
    field = variable === :temperature ? hydro.temperature : hydro.salinity
    ts = Float64.(field[i_st, j_st, k_surf])
    times_hr = hydro.time_seconds / 3600.0

    fig = Figure(size = (600, 360), fontsize = 12)
    ax = Axis(fig[1, 1],
        title = "Eulerian $(variable) at Station\n$(round(lons[i_st], digits=2))°E, $(round(lats[j_st], digits=2))°N",
        xlabel = "Time (h)",
        ylabel = variable === :temperature ? "Temperature (°C)" : "Salinity (PSU)"
    )
    lines!(ax, [times_hr], [ts]; color = variable === :temperature ? :crimson : :navy, linewidth = 2)
    scatter!(ax, [times_hr], [ts]; color = variable === :temperature ? :crimson : :navy, markersize = 4)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_hovmoller(
        hydro::NamedTuple;
        station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        variable::Symbol = :temperature,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_hovmoller.png"
    ) -> Figure

Plot Hovmöller diagram (depth vs time) at a station.
"""
function plot_hydrodynamic_hovmoller(
    hydro::NamedTuple;
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    variable::Symbol = :temperature,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_hovmoller.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths

    i_st = isnothing(station) ? length(lons) ÷ 2 : argmin(abs.(lons .- station[1]))
    j_st = isnothing(station) ? length(lats) ÷ 2 : argmin(abs.(lats .- station[2]))

    field = variable === :temperature ? hydro.temperature : hydro.salinity
    nz = length(depths)
    n_t = hydro isa NamedTuple && hasproperty(hydro, :times) ? length(hydro.times) : 1

    # For single-interval data, create a 1-column Hovmöller
    sec = reshape(Float64.(field[i_st, j_st, :]), 1, nz)

    fig = Figure(size = (420, 480), fontsize = 12)
    ax = Axis(fig[1, 1],
        title = "$(variable) Hovmöller at Station",
        xlabel = "Time",
        ylabel = "Depth (m)"
    )
    hm = heatmap!(ax, [1], depths, sec;
                  colormap = :thermal,
                  colorrange = _finite_colorrange(sec; fallback = (0.0, 15.0)))
    Colorbar(fig[1, 2], hm, label = variable === :temperature ? "°C" : "PSU")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_bottom_temperature(
        hydro::NamedTuple;
        output_path::Union{Nothing, AbstractString} = "outputs/bottom_temperature_map.png"
    ) -> Figure

Plot bottom temperature map.
"""
function plot_hydrodynamic_bottom_temperature(
    hydro::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/bottom_temperature_map.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths

    # Near-bed temperature: deepest finite value in each water column
    nx, ny, nz = size(hydro.temperature)
    order_deep_first = sortperm(depths)
    T_bottom = fill(NaN, nx, ny)
    for j in 1:ny, i in 1:nx
        for k in order_deep_first
            val = hydro.temperature[i, j, k]
            if isfinite(val)
                T_bottom[i, j] = Float64(val)
                break
            end
        end
    end

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Bottom Temperature (°C)",
        xlabel = "Longitude (°E)",
        ylabel = "Latitude (°N)",
        aspect = DataAspect()
    )
    hm = heatmap!(ax, lons, lats, T_bottom, colormap = :thermal;
                   colorrange = _finite_colorrange(T_bottom; fallback = (0.0, 15.0)))
    Colorbar(fig[1, 2], hm, label = "°C")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_ts_diagram(
        hydro::NamedTuple;
        output_path::Union{Nothing, AbstractString} = "outputs/ts_diagram.png"
    ) -> Figure

Plot T-S diagram (temperature vs salinity scatter).
"""
function plot_hydrodynamic_ts_diagram(
    hydro::NamedTuple;
    output_path::Union{Nothing, AbstractString} = "outputs/ts_diagram.png"
)
    T_all = Float64.(hydro.temperature[:])
    S_all = Float64.(hydro.salinity[:])

    # Filter non-finite values
    valid = isfinite.(T_all) .& isfinite.(S_all)

    fig = Figure(size = (500, 450), fontsize = 12)
    ax = Axis(fig[1, 1],
        title = "T-S Diagram",
        xlabel = "Salinity (PSU)",
        ylabel = "Temperature (°C)"
    )
    scatter!(ax, S_all[valid], T_all[valid]; color = :steelblue, markersize = 3, alpha = 0.5)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end



