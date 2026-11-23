"""
    larval_disperal.jl

Larval dispersal figures.

Trajectory maps, vertical-migration profiles, run comparison, settlement
density, empirical movement fields, connectivity and recruitment summaries.
"""

using CairoMakie
import CairoMakie: record

"""
    plot_particle_trajectories(
        trajectories::NamedTuple;
        bathymetry_data::Union{Nothing, NamedTuple} = nothing,
        strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
        title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
        output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
        max_display_particles::Int = 100
    )

Render a 2D spatial map of snow crab larval drift trajectories with optional
background bathymetry and management strata boundaries.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort` with fields
  `(lons, lats, depths, times, ids)`.
- `bathymetry_data::Union{Nothing, NamedTuple}`: Optional `(lon, lat, elevation)`
  for background contours.
- `strata::Union{Nothing, AbstractVector{<:NamedTuple}}`: Optional management strata
  definitions with polygon coordinates.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save the figure (PNG/PDF).
- `max_display_particles::Int`: Maximum number of individual tracks to render.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- North, E. W., et al. (2009). ICES Cooperative Research Report, No. 295.
"""
function plot_particle_trajectories(
    trajectories::NamedTuple;
    bathymetry_data::Union{Nothing, NamedTuple} = nothing,
    strata::Union{Nothing, AbstractVector{<:NamedTuple}} = nothing,
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::AbstractString = "Snow Crab Larval Dispersal Trajectories",
    output_path::Union{Nothing, AbstractString} = "outputs/particle_trajectories.png",
    max_display_particles::Int = 100
)
    use_rot = abs(Float64(rotation_deg)) > 1e-3
    bathy_mp = nothing
    shelf_frame = nothing
    if !isnothing(bathymetry_data)
        elev = Float64.(bathymetry_data.elevation)
        elev[elev .>= 0.0] .= NaN
        bathy_mp = _prepare_map_fields(
            bathymetry_data.lon, bathymetry_data.lat, elev, elev;
            rotation_deg = rotation_deg
        )
        shelf_frame = bathy_mp.shelf_frame
    elseif use_rot
        fin_lo = filter(isfinite, vec(trajectories.lons))
        fin_la = filter(isfinite, vec(trajectories.lats))
        lon0 = isempty(fin_lo) ? -62.5 : mean(fin_lo)
        lat0 = isempty(fin_la) ? 44.2 : mean(fin_la)
        shelf_frame = (Float64(lon0), Float64(lat0), Float64(rotation_deg))
    end

    proj_pt(lo, la) = isnothing(shelf_frame) ?
        (Float64(lo), Float64(la)) :
        _lonlat_to_shelf_km(lo, la, shelf_frame[1], shelf_frame[2], shelf_frame[3])

    xlabel_txt = isnothing(shelf_frame) ? "Longitude (°E)" :
        "Along-Shelf Distance (km, $(round(Int, shelf_frame[3]))° CCW from E)"
    ylabel_txt = isnothing(shelf_frame) ? "Latitude (°N)" : "Cross-Shelf Distance (km)"

    fig = Figure(size = (980, 620), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
        aspect = DataAspect()
    )

    # 1. Plot background bathymetry if provided (mask land >= 0 to NaN so land stays background)
    if !isnothing(bathy_mp) && any(isfinite, bathy_mp.fields[1])
        co = contourf!(
            ax,
            bathy_mp.xs,
            bathy_mp.ys,
            bathy_mp.fields[1],
            colormap = :viridis,
            levels = 15
        )
        Colorbar(fig[1, 2], co, label = "Seafloor Elevation (m)")
    end

    # 2. Draw management strata outlines if provided
    if !isnothing(strata)
        for s in strata
            poly_lons = Float64[]
            poly_lats = Float64[]
            if hasproperty(s, :lons) && hasproperty(s, :lats) && !isempty(s.lons)
                poly_lons = Float64.(s.lons)
                poly_lats = Float64.(s.lats)
            elseif hasproperty(s, :polygon) && !isempty(s.polygon)
                poly_lons = [Float64(pt[1]) for pt in s.polygon]
                poly_lats = [Float64(pt[2]) for pt in s.polygon]
            end
            if !isempty(poly_lons)
                # Close the ring: an unclosed polygon renders as a gap, and a gap in
                # a management boundary reads as an opening that is not there.
                if poly_lons[1] != poly_lons[end] || poly_lats[1] != poly_lats[end]
                    push!(poly_lons, poly_lons[1])
                    push!(poly_lats, poly_lats[1])
                end
                pts_s = [proj_pt(poly_lons[i], poly_lats[i]) for i in eachindex(poly_lons)]
                lines!(
                    ax,
                    first.(pts_s),
                    last.(pts_s),
                    color = (:gold, 0.75),
                    linewidth = 1.8,
                    linestyle = :dash
                )
            end
        end
    end

    # 3. Plot larval tracks with high-visibility outline underlay and vivid color
    n_particles = size(trajectories.lons, 1)
    n_plot = min(n_particles, max_display_particles)
    track_xs = Vector{Vector{Float64}}(undef, n_plot)
    track_ys = Vector{Vector{Float64}}(undef, n_plot)
    for p in 1:n_plot
        nt_p = size(trajectories.lons, 2)
        xs_p = Vector{Float64}(undef, nt_p)
        ys_p = Vector{Float64}(undef, nt_p)
        for t in 1:nt_p
            lo, la = trajectories.lons[p, t], trajectories.lats[p, t]
            if isfinite(lo) && isfinite(la)
                xs_p[t], ys_p[t] = proj_pt(lo, la)
            else
                xs_p[t], ys_p[t] = NaN, NaN
            end
        end
        track_xs[p] = xs_p
        track_ys[p] = ys_p
    end

    # High-contrast dark shadow underlay
    for p in 1:n_plot
        lines!(ax, track_xs[p], track_ys[p], color = (:black, 0.75), linewidth = 3.8)
    end

    # Vivid foreground trajectory lines
    for p in 1:n_plot
        lines!(ax, track_xs[p], track_ys[p], color = :cyan, linewidth = 2.2)
    end

    # 4. Mark release positions (start) and final drift positions (end) with strokes
    rel_x = [track_xs[p][1] for p in 1:n_plot]
    rel_y = [track_ys[p][1] for p in 1:n_plot]
    end_x = [track_xs[p][end] for p in 1:n_plot]
    end_y = [track_ys[p][end] for p in 1:n_plot]
    scatter!(
        ax, rel_x, rel_y,
        color = :springgreen,
        strokecolor = :black,
        strokewidth = 1.5,
        marker = :circle,
        markersize = 11,
        label = "Release (t = 0)"
    )
    scatter!(
        ax, end_x, end_y,
        color = :crimson,
        strokecolor = :white,
        strokewidth = 1.5,
        marker = :diamond,
        markersize = 12,
        label = "Final Position"
    )

    _overlay_coastline!(ax, coastline; shelf_frame = shelf_frame)
    if !isnothing(bathy_mp)
        xlims!(ax, bathy_mp.xs[1], bathy_mp.xs[end])
        ylims!(ax, bathy_mp.ys[1], bathy_mp.ys[end])
    elseif !isnothing(coastline)
        fin_x = filter(isfinite, reduce(vcat, track_xs; init = Float64[]))
        fin_y = filter(isfinite, reduce(vcat, track_ys; init = Float64[]))
        if !isempty(fin_x) && !isempty(fin_y)
            pad = isnothing(shelf_frame) ? 0.5 : 25.0
            xlims!(ax, minimum(fin_x) - pad, maximum(fin_x) + pad)
            ylims!(ax, minimum(fin_y) - pad, maximum(fin_y) + pad)
        end
    end

    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_vertical_migration_profiles(
        trajectories::NamedTuple;
        sample_indices = 1:min(5, size(trajectories.depths, 1)),
        title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
        show_surface_target::Bool = true,
        surface_target_depth::Real = -10.0,
        output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
    )

Plot larval depth trajectories over time to visualize initial ascent from benthic
release depths and diurnal oscillations between nighttime surface grazing and
daytime sub-surface predator avoidance.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `sample_indices`: Range or vector of particle indices to plot.
- `title::AbstractString`: Figure title.
- `show_surface_target::Bool`: Whether to draw reference line for surface target layer.
- `surface_target_depth::Real`: Surface mixed layer target depth in meters (default -10.0 m).
- `output_path::Union{Nothing, AbstractString}`: Destination path for figure.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Incze, L. S., et al. (1987). *Marine Biology*, 95(2), 195-200.
- Sainte-Marie, G., & Sainte-Marie, B. (1999). *Can. J. Fish. Aquat. Sci.*, 56(11), 2181-2193.
"""
function plot_vertical_migration_profiles(
    trajectories::NamedTuple;
    sample_indices = 1:min(5, size(trajectories.depths, 1)),
    title::AbstractString = "Diel Vertical Migration (DVM) Depth Profiles",
    show_surface_target::Bool = true,
    surface_target_depth::Real = -10.0,
    output_path::Union{Nothing, AbstractString} = "outputs/dvm_profiles.png"
)
    fig = Figure(size = (850, 480), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Simulation Time (hours)",
        ylabel = "Depth (m)"
    )

    time_hours = trajectories.times ./ 3600.0
    colors = [:royalblue, :darkorange, :forestgreen, :purple, :crimson]

    for (idx, p) in enumerate(sample_indices)
        c = colors[mod1(idx, length(colors))]
        lines!(
            ax,
            time_hours,
            trajectories.depths[p, :],
            color = c,
            linewidth = 2.0,
            label = "Larva $(p)"
        )
        scatter!(
            ax,
            [time_hours[1]],
            [trajectories.depths[p, 1]],
            color = c,
            marker = :circle,
            markersize = 8
        )
    end

    hlines!(ax, [0.0], color = :gray50, linestyle = :dash)
    if show_surface_target
        hlines!(
            ax,
            [Float64(surface_target_depth)],
            color = :teal,
            linestyle = :dot,
            linewidth = 1.5,
            label = "Surface Target ($(round(surface_target_depth, digits=1))m)"
        )
    end
    axislegend(ax, position = :rb)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    compare_run_dispersal(
        run_dict::Union{AbstractDict, NamedTuple};
        title::AbstractString = "Larval Dispersal Across Runs",
        output_path::Union{Nothing, AbstractString} = "outputs/run_comparison.png"
    )

Compare spatial trajectories and final dispersion centroids across several
archived runs.

# Mathematical Formulations
- **Centroid Center of Mass**:
  ```math
  \\bar{\\\\\\\\lambda} = \\\\\\\\\\\\\\frac{1}{N_p} \\sum_{p=1}^{N_p} \\\\\\\\lambda_p(t_{\\text{end}}), \\quad
  \\bar{\\phi} = \\\\\\\\\\\\\\frac{1}{N_p} \\sum_{p=1}^{N_p} \\phi_p(t_{\\text{end}})
  ```

# Inputs
- `run_dict::Union{AbstractDict, NamedTuple}`: Mapping of run identifiers or keys
  to trajectory NamedTuples (or data containers holding a `.trajectories` field).
- `title::AbstractString`: Comparison figure title.
- `output_path::Union{Nothing, AbstractString}`: Destination filepath to save output.

# Outputs
- `Figure`: CairoMakie figure object.

# References
- Brickman, D., et al. (2018). *Progress in Oceanography*, 164, 49-64.
"""
function compare_run_dispersal(
    run_dict::Union{AbstractDict, NamedTuple};
    title::AbstractString = "Larval Dispersal Across Runs",
    output_path::Union{Nothing, AbstractString} = "outputs/run_comparison.png"
)
    n_runs = length(run_dict)
    n_cols = min(n_runs, 3)
    n_rows = ceil(Int, n_runs / 3)
    fig = Figure(size = (450 * n_cols, 420 * n_rows), fontsize = 12)

    run_keys = run_dict isa NamedTuple ? keys(run_dict) : collect(keys(run_dict))
    palette = (:royalblue, :seagreen, :darkorange, :orangered, :crimson, :firebrick)

    for (r_idx, r_key) in enumerate(run_keys)
        row = div(r_idx - 1, 3) + 1
        col = mod(r_idx - 1, 3) + 1

        raw_val = run_dict isa NamedTuple ? getproperty(run_dict, r_key) : run_dict[r_key]
        traj = hasproperty(raw_val, :trajectories) ? raw_val.trajectories : raw_val
        c = palette[mod1(r_idx, length(palette))]

        ax = Axis(
            fig[row, col],
            title = "$(string(r_key))",
            xlabel = "Longitude (°E)",
            ylabel = "Latitude (°N)"
        )

        n_p = size(traj.lons, 1)
        for p in 1:min(n_p, 50)
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (:black, 0.5), linewidth = 2.8
            )
            lines!(
                ax, traj.lons[p, :], traj.lats[p, :],
                color = (c, 0.85), linewidth = 1.8
            )
        end

        mean_lon_end = sum(traj.lons[:, end]) / n_p
        mean_lat_end = sum(traj.lats[:, end]) / n_p

        scatter!(
            ax,
            [mean_lon_end],
            [mean_lat_end],
            color = :gold,
            strokecolor = :black,
            strokewidth = 1.5,
            marker = :star5,
            markersize = 16,
            label = "Centroid"
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_larval_dispersal_density(
        trajectories::NamedTuple;
        lon_bins = range(-68.0, -57.0, length = 100),
        lat_bins = range(41.0, 48.0, length = 100),
        title::AbstractString = "Larval Settlement Density Kernel",
        output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
    )

Compute and plot a 2D histogram heatmap of final particle positions representing
larval settlement and nursery retention hotspots.

# Inputs
- `trajectories::NamedTuple`: Output from `track_larval_cohort`.
- `lon_bins`: Grid longitude bin edges.
- `lat_bins`: Grid latitude bin edges.
- `title::AbstractString`: Figure title.
- `output_path::Union{Nothing, AbstractString}`: Path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_larval_dispersal_density(
    trajectories::NamedTuple;
    lon_bins = range(-68.0, -57.0, length = 100),
    lat_bins = range(41.0, 48.0, length = 100),
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::AbstractString = "Larval Settlement Density Kernel",
    output_path::Union{Nothing, AbstractString} = "outputs/settlement_density.png"
)
    end_lons = trajectories.lons[:, end]
    end_lats = trajectories.lats[:, end]
    n_p = length(end_lons)

    nx = length(lon_bins) - 1
    ny = length(lat_bins) - 1
    density = zeros(Float64, nx, ny)

    for p in 1:n_p
        x = end_lons[p]
        y = end_lats[p]
        if lon_bins[1] <= x <= lon_bins[end] && lat_bins[1] <= y <= lat_bins[end]
            i = clamp(searchsortedlast(lon_bins, x), 1, nx)
            j = clamp(searchsortedlast(lat_bins, y), 1, ny)
            density[i, j] += 1.0
        end
    end

    density ./= max(1.0, sum(density))
    density .*= 100.0
    density[density .<= 0.0] .= NaN

    lon_centers = [(lon_bins[i] + lon_bins[i + 1]) / 2.0 for i in 1:nx]
    lat_centers = [(lat_bins[j] + lat_bins[j + 1]) / 2.0 for j in 1:ny]
    mp = _prepare_map_fields(
        lon_centers, lat_centers, density, density;
        rotation_deg = rotation_deg, nearest = true
    )
    dens_p = mp.fields[1]

    fig = Figure(size = (950, 560), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm = heatmap!(
        ax, mp.xs, mp.ys, dens_p;
        colormap = :inferno,
        colorrange = _finite_colorrange(dens_p; fallback = (0.0, 1.0)),
        nan_color = :transparent
    )
    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax, mp.xs[1], mp.xs[end])
    ylims!(ax, mp.ys[1], mp.ys[end])
    Colorbar(fig[1, 2], hm, label = "Settlement Density (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_empirical_movement_field(
        emp_mov::NamedTuple;
        title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
        output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
    )

Plot empirical velocity vector arrows over the empirical turbulent diffusivity field.

# Inputs
- `emp_mov::NamedTuple`: Output from `estimate_empirical_movement`.
- `title::AbstractString`: Plot title.
- `output_path::Union{Nothing, AbstractString}`: File path to save figure.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_empirical_movement_field(
    emp_mov::NamedTuple;
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::AbstractString = "Empirical Advection & Turbulent Diffusivity Field",
    output_path::Union{Nothing, AbstractString} = "outputs/empirical_movement.png"
)
    lon_c = emp_mov.lon_centers
    lat_c = emp_mov.lat_centers
    diff = Float64.(emp_mov.diffusivity)
    diff[(emp_mov.sample_count .<= 0) .| (diff .<= 0.0)] .= NaN

    mp = _prepare_map_fields(
        lon_c, lat_c, diff, diff;
        rotation_deg = rotation_deg, nearest = true
    )
    diff_p = mp.fields[1]

    fig = Figure(size = (950, 560), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )

    hm = heatmap!(
        ax, mp.xs, mp.ys, diff_p;
        colormap = :viridis,
        colorrange = _finite_colorrange(diff_p; fallback = (0.0, 1.0)),
        nan_color = :transparent
    )
    Colorbar(fig[1, 2], hm, label = "Empirical Diffusivity D (m² s⁻¹)")

    pts = Point2f[]
    dirs = Vec2f[]
    u_m = emp_mov.u_mean
    v_m = emp_mov.v_mean

    for i in 1:length(lon_c), j in 1:length(lat_c)
        if !isnan(u_m[i, j]) && !isnan(v_m[i, j]) && emp_mov.sample_count[i, j] > 2
            if isnothing(mp.shelf_frame)
                push!(pts, Point2f(lon_c[i], lat_c[j]))
                push!(dirs, Vec2f(u_m[i, j] * 0.5, v_m[i, j] * 0.5))
            else
                lon0, lat0, rdeg = mp.shelf_frame
                xk, yk = _lonlat_to_shelf_km(lon_c[i], lat_c[j], lon0, lat0, rdeg)
                c_r, s_r = cosd(rdeg), sind(rdeg)
                ur = u_m[i, j] * c_r + v_m[i, j] * s_r
                vr = -u_m[i, j] * s_r + v_m[i, j] * c_r
                push!(pts, Point2f(xk, yk))
                push!(dirs, Vec2f(ur * 45.0, vr * 45.0))
            end
        end
    end

    if !isempty(pts)
        if isdefined(Makie, :arrows2d!)
            arrows2d!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        else
            arrows!(
                ax,
                pts,
                dirs,
                tipwidth = 6,
                tiplength = 10,
                lengthscale = 1.0,
                tipcolor = :white,
                shaftcolor = :white
            )
        end
    end

    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax, mp.xs[1], mp.xs[end])
    ylims!(ax, mp.ys[1], mp.ys[end])

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_connectivity_matrix(
        conn::NamedTuple;
        title::AbstractString = "Macro-Regional Population Connectivity Matrix",
        output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
    )

Render an annotated transition probability heatmap showing larval transfer probabilities
and self-retention percentages between management zones.

# Inputs
- `conn::NamedTuple`: Output from `compute_empirical_connectivity`.
- `title::AbstractString`: Title string.
- `output_path::Union{Nothing, AbstractString}`: Path for output image.

# Outputs
- `Figure`: CairoMakie figure object.
"""
function plot_connectivity_matrix(
    conn::NamedTuple;
    title::AbstractString = "Macro-Regional Population Connectivity Matrix",
    output_path::Union{Nothing, AbstractString} = "outputs/connectivity_matrix.png"
)
    mat = conn.matrix .* 100.0
    strata = conn.strata_names
    n_s = length(strata)

    fig = Figure(size = (800, 650), fontsize = 12)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Destination Stratum",
        ylabel = "Source Stratum",
        xticks = (1:n_s, strata),
        yticks = (1:n_s, strata),
        xticklabelrotation = 0.35
    )

    hm = heatmap!(ax, 1:n_s, 1:n_s, mat', colormap = :blues)
    Colorbar(fig[1, 2], hm, label = "Transition Probability (%)")

    for i in 1:n_s, j in 1:n_s
        val = mat[i, j]
        txt_col = val > 50.0 ? :white : :black
        text!(
            ax,
            Point2f(j, i),
            text = string(round(val, digits = 1), "%"),
            align = (:center, :center),
            color = txt_col,
            fontsize = 12
        )
    end

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_thermal_exposure_map(
        therm_met::NamedTuple;
        title::AbstractString = "Larval Thermal Exposure & Degree-Days",
        output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
    )

Plot a two-panel map showing accumulated thermal degree-days and mean exposure temperature.
"""
function plot_thermal_exposure_map(
    therm_met::NamedTuple;
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::AbstractString = "Larval Thermal Exposure & Degree-Days",
    output_path::Union{Nothing, AbstractString} = "outputs/thermal_exposure.png"
)
    lon_c = therm_met.lon_centers
    lat_c = therm_met.lat_centers
    dd = Float64.(therm_met.mean_degree_days)
    temp = Float64.(therm_met.mean_exposure_temperature)
    mp = _prepare_map_fields(
        lon_c, lat_c, dd, dd, temp;
        rotation_deg = rotation_deg, nearest = true
    )
    dd_p, temp_p = mp.fields[1], mp.fields[2]

    fig = Figure(size = (1200, 480), fontsize = 12)

    # Panel 1: Accumulated Degree-Days
    ax1 = Axis(
        fig[1, 1],
        title = "Mean Cumulative Degree-Days (°C·days)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm1 = heatmap!(
        ax1, mp.xs, mp.ys, dd_p;
        colormap = :inferno,
        colorrange = _finite_colorrange(dd_p; fallback = (0.0, 1.0)),
        nan_color = :transparent
    )
    _overlay_coastline!(ax1, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax1, mp.xs[1], mp.xs[end])
    ylims!(ax1, mp.ys[1], mp.ys[end])
    Colorbar(fig[1, 2], hm1, label = "Degree-Days (°C·d)")

    # Panel 2: Mean Ambient Exposure Temperature
    ax2 = Axis(
        fig[1, 3],
        title = "Mean Exposure Temperature (°C)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm2 = heatmap!(
        ax2, mp.xs, mp.ys, temp_p;
        colormap = :thermal,
        colorrange = _finite_colorrange(temp_p; fallback = (0.0, 15.0)),
        nan_color = :transparent
    )
    _overlay_coastline!(ax2, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax2, mp.xs[1], mp.xs[end])
    ylims!(ax2, mp.ys[1], mp.ys[end])
    Colorbar(fig[1, 4], hm2, label = "Temperature (°C)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_recruitment_summary(
        rec_met::NamedTuple;
        title::AbstractString = "Larval Recruitment & Retention Summary",
        output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
    )

Plot a multi-panel spatial distribution of spawning release, settlement density, and recruitment success rate.
"""
function plot_recruitment_summary(
    rec_met::NamedTuple;
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::AbstractString = "Larval Recruitment & Retention Summary",
    output_path::Union{Nothing, AbstractString} = "outputs/recruitment_summary.png"
)
    lon_c = rec_met.lon_centers
    lat_c = rec_met.lat_centers
    dens = Float64.(rec_met.settlement_density)
    dens[dens .<= 0.0] .= NaN
    succ_pct = Float64.(rec_met.success_rate) .* 100.0
    succ_pct[(rec_met.settlement_density .<= 0.0) .& (succ_pct .<= 0.0)] .= NaN

    mp = _prepare_map_fields(
        lon_c, lat_c, dens, dens, succ_pct;
        rotation_deg = rotation_deg, nearest = true
    )
    dens_p, succ_p = mp.fields[1], mp.fields[2]

    fig = Figure(size = (1200, 480), fontsize = 12)

    # Panel 1: Settlement Density
    ax1 = Axis(
        fig[1, 1],
        title = "Settlement Density (count)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm1 = heatmap!(
        ax1, mp.xs, mp.ys, dens_p;
        colormap = :plasma,
        colorrange = _finite_colorrange(dens_p; fallback = (0.0, 1.0)),
        nan_color = :transparent
    )
    _overlay_coastline!(ax1, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax1, mp.xs[1], mp.xs[end])
    ylims!(ax1, mp.ys[1], mp.ys[end])
    Colorbar(fig[1, 2], hm1, label = "Settled Larvae")

    # Panel 2: Recruitment Success Rate
    ax2 = Axis(
        fig[1, 3],
        title = "Benthic Recruitment Success Rate (%)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm2 = heatmap!(
        ax2, mp.xs, mp.ys, succ_p;
        colormap = :viridis,
        colorrange = _finite_colorrange(succ_p; fallback = (0.0, 100.0)),
        nan_color = :transparent
    )
    _overlay_coastline!(ax2, coastline; shelf_frame = mp.shelf_frame)
    xlims!(ax2, mp.xs[1], mp.xs[end])
    ylims!(ax2, mp.ys[1], mp.ys[end])
    Colorbar(fig[1, 4], hm2, label = "Success Rate (%)")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_bottom_temperature_map(hydrodynamics; depth = nothing, output_path = ...)

Map of near-bed temperature across the domain.

# Inputs
- `hydrodynamics`: Extracted hydrodynamic `NamedTuple` or Zarr store path.
- `depth`: Optional target depth (m); defaults to the deepest resolved level.
- `domain_lon`, `domain_lat`: Longitude and latitude bounds.
- `title`: Optional figure title.
- `output_path`: PNG path, or `nothing`.

# Outputs
- `Makie.Figure`: Near-bed temperature map.
"""
function plot_bottom_temperature_map(
    hydrodynamics::Any;
    depth::Union{Nothing, Real} = nothing,
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    coastline = nothing,
    rotation_deg::Real = 0.0,
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/bottom_temperature_map.png"
)
    f = extract_hydrodynamic_dataset(
        hydrodynamics;
        domain_lon = domain_lon,
        domain_lat = domain_lat
    )
    zs = collect(Float64, f.depths)
    k = isnothing(depth) ? argmin(zs) : argmin(abs.(zs .- Float64(depth)))
    Tbot = Float64.(f.temperature[:, :, k])
    if hasproperty(f, :bathymetry) && !isnothing(f.bathymetry)
        mask_bathy!(Tbot, f.bathymetry, k)
    end
    mp = _prepare_map_fields(f.lons, f.lats, Tbot, Tbot; rotation_deg = rotation_deg)
    Tbot_p = mp.fields[1]
    cr = _finite_colorrange(Tbot_p; fallback = (0.0, 1.0))
    fig = Figure(size = (950, 520), fontsize = 12)
    ax = Axis(
        fig[1, 1],
        title = something(
            title,
            "Near-Bed Temperature (deepest resolved level, $(zs[k]) m)"
        ),
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm = heatmap!(
        ax, mp.xs, mp.ys, Tbot_p;
        colormap = :thermal,
        colorrange = cr,
        nan_color = :transparent
    )
    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax, mp.xs, mp.ys, Tbot_p)
    else
        xlims!(ax, mp.xs[1], mp.xs[end]); ylims!(ax, mp.ys[1], mp.ys[end])
    end
    Colorbar(fig[1, 2], hm, label = "Temperature (°C)")
    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

"""
    plot_temperature_salinity_diagram(hydrodynamics; output_path = ...)

Temperature–salinity (T–S) scatter of the domain, coloured by depth.

# Inputs
- `hydrodynamics`: Extracted hydrodynamic `NamedTuple` or Zarr store path.
- `domain_lon`, `domain_lat`: Longitude and latitude bounds.
- `title`: Optional figure title.
- `output_path`: PNG path, or `nothing`.

# Outputs
- `Makie.Figure`: Temperature–salinity scatter diagram.
"""
function plot_temperature_salinity_diagram(
    hydrodynamics::Any;
    domain_lon::Tuple{<:Real, <:Real} = (-71.0, -53.0),
    domain_lat::Tuple{<:Real, <:Real} = (40.0, 48.5),
    title::Union{Nothing, AbstractString} = nothing,
    output_path::Union{Nothing, AbstractString} = "outputs/ts_diagram.png"
)
    f = extract_hydrodynamic_dataset(
        hydrodynamics;
        domain_lon = domain_lon,
        domain_lat = domain_lat
    )
    zs = collect(Float64, f.depths)
    nx, ny, _ = size(f.salinity)
    Sv = vec(f.salinity)
    Tv = vec(f.temperature)
    Zv = repeat(zs, inner = nx * ny)

    fig = Figure(size = (760, 560), fontsize = 12)
    ax = Axis(
        fig[1, 1],
        title = something(title, "Temperature–Salinity Diagram"),
        xlabel = "Salinity (psu)",
        ylabel = "Temperature (°C)"
    )
    m = isfinite.(Tv) .& isfinite.(Sv)
    sc = scatter!(ax, Sv[m], Tv[m]; color = Zv[m], colormap = :viridis)
    Colorbar(fig[1, 2], sc, label = "Depth (m)")
    isnothing(output_path) || (mkpath(dirname(output_path)); save(output_path, fig))
    return fig
end

