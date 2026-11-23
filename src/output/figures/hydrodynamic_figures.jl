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
    _lonlat_to_shelf_km(lon, lat, lon0, lat0, rot_deg) -> Tuple{Float64, Float64}

Project `(lon, lat)` into a local tangent-plane frame centered at `(lon0, lat0)` rotated
`rot_deg` degrees counter-clockwise from East (along-shelf x, cross-shelf y in km).
"""
@inline function _lonlat_to_shelf_km(lon::Real, lat::Real, lon0::Real, lat0::Real, rot_deg::Real)
    km_per_lat = 111.132
    km_per_lon = 111.320 * cosd(Float64(lat0))
    dx = (Float64(lon) - Float64(lon0)) * km_per_lon
    dy = (Float64(lat) - Float64(lat0)) * km_per_lat
    c = cosd(Float64(rot_deg))
    s = sind(Float64(rot_deg))
    return (dx * c + dy * s, -dx * s + dy * c)
end

@inline function _shelf_km_to_lonlat(x_km::Real, y_km::Real, lon0::Real, lat0::Real, rot_deg::Real)
    km_per_lat = 111.132
    km_per_lon = 111.320 * cosd(Float64(lat0))
    c = cosd(Float64(rot_deg))
    s = sind(Float64(rot_deg))
    dx = Float64(x_km) * c - Float64(y_km) * s
    dy = Float64(x_km) * s + Float64(y_km) * c
    return (Float64(lon0) + dx / km_per_lon, Float64(lat0) + dy / km_per_lat)
end

function _bilinear_2d(lons::AbstractVector, lats::AbstractVector, A::AbstractMatrix, lo::Float64, la::Float64)::Float64
    nx, ny = size(A)
    (lo < lons[1] || lo > lons[end] || la < lats[1] || la > lats[end]) && return NaN
    i1 = clamp(searchsortedfirst(lons, lo), 2, nx)
    j1 = clamp(searchsortedfirst(lats, la), 2, ny)
    i0 = i1 - 1
    j0 = j1 - 1
    sx = lons[i1] == lons[i0] ? 0.0 : clamp((lo - lons[i0]) / (lons[i1] - lons[i0]), 0.0, 1.0)
    sy = lats[j1] == lats[j0] ? 0.0 : clamp((la - lats[j0]) / (lats[j1] - lats[j0]), 0.0, 1.0)
    @inbounds begin
        v00 = Float64(A[i0, j0]); v10 = Float64(A[i1, j0])
        v01 = Float64(A[i0, j1]); v11 = Float64(A[i1, j1])
    end
    res = (1 - sx) * (1 - sy) * v00 + sx * (1 - sy) * v10 +
          (1 - sx) * sy * v01 + sx * sy * v11
    return res
end

function _nearest_2d(lons::AbstractVector, lats::AbstractVector, A::AbstractMatrix, lo::Float64, la::Float64)::Float64
    nx, ny = size(A)
    (lo < lons[1] || lo > lons[end] || la < lats[1] || la > lats[end]) && return NaN
    i1 = clamp(searchsortedfirst(lons, lo), 1, nx)
    i = (i1 > 1 && abs(lo - lons[i1 - 1]) < abs(lons[i1] - lo)) ? (i1 - 1) : i1
    j1 = clamp(searchsortedfirst(lats, la), 1, ny)
    j = (j1 > 1 && abs(la - lats[j1 - 1]) < abs(lats[j1] - la)) ? (j1 - 1) : j1
    return @inbounds Float64(A[i, j])
end

"""
    _resample_to_shelf_grid(lons, lats, ref_field, fields...; rot_deg = 35.0, n_along = 260, n_cross = 180, nearest = false)

Resample one or more 2D `(lon, lat)` fields onto a regular shelf-aligned `(x_km, y_km)`
grid tightly bounding the finite wet cells, eliminating empty rotated corner whitespace.
Use `nearest = true` for sparse binned larval grids so isolated valid cells are preserved.
"""
function _resample_to_shelf_grid(
    lons::AbstractVector,
    lats::AbstractVector,
    ref_field::AbstractMatrix,
    fields::AbstractMatrix...;
    rot_deg::Real = 35.0,
    n_along::Int = 260,
    n_cross::Int = 180,
    nearest::Bool = false
)
    wet_ij = findall(isfinite, ref_field)
    lon0 = isempty(wet_ij) ? 0.5 * (lons[1] + lons[end]) : mean(lons[idx[1]] for idx in wet_ij)
    lat0 = isempty(wet_ij) ? 0.5 * (lats[1] + lats[end]) : mean(lats[idx[2]] for idx in wet_ij)

    x_min, x_max = Inf, -Inf
    y_min, y_max = Inf, -Inf
    for idx in wet_ij
        xk, yk = _lonlat_to_shelf_km(lons[idx[1]], lats[idx[2]], lon0, lat0, rot_deg)
        xk < x_min && (x_min = xk); xk > x_max && (x_max = xk)
        yk < y_min && (y_min = yk); yk > y_max && (y_max = yk)
    end
    if !isfinite(x_min)
        x_min, x_max = -200.0, 200.0
        y_min, y_max = -150.0, 150.0
    end
    xs = collect(range(x_min - 15.0, x_max + 15.0; length = n_along))
    ys = collect(range(y_min - 15.0, y_max + 15.0; length = n_cross))
    outs = [fill(NaN, n_along, n_cross) for _ in 1:length(fields)]
    lons_v = Float64.(collect(lons))
    lats_v = Float64.(collect(lats))
    for j in 1:n_cross, i in 1:n_along
        lo, la = _shelf_km_to_lonlat(xs[i], ys[j], lon0, lat0, rot_deg)
        for (m, F) in enumerate(fields)
            outs[m][i, j] = nearest ?
                _nearest_2d(lons_v, lats_v, F, lo, la) :
                _bilinear_2d(lons_v, lats_v, F, lo, la)
        end
    end
    return (xs = xs, ys = ys, fields = outs, lon0 = lon0, lat0 = lat0, rot_deg = Float64(rot_deg))
end

"""
    _prepare_map_fields(lons, lats, ref_field, fields...; rotation_deg = 0.0, nearest = false)

Return `(xs, ys, plot_fields, xlabel, ylabel, shelf_frame)` in either native `(lon, lat)`
coordinates (`rotation_deg == 0`) or the rotated shelf-aligned km frame.
"""
function _prepare_map_fields(
    lons::AbstractVector,
    lats::AbstractVector,
    ref_field::AbstractMatrix,
    fields::AbstractMatrix...;
    rotation_deg::Real = 0.0,
    nearest::Bool = false
)
    if abs(Float64(rotation_deg)) > 1e-3
        res = _resample_to_shelf_grid(
            lons, lats, ref_field, fields...;
            rot_deg = rotation_deg, nearest = nearest
        )
        return (
            xs = res.xs,
            ys = res.ys,
            fields = res.fields,
            xlabel = "Along-Shelf Distance (km, $(round(Int, res.rot_deg))° CCW from E)",
            ylabel = "Cross-Shelf Distance (km)",
            shelf_frame = (res.lon0, res.lat0, res.rot_deg)
        )
    else
        return (
            xs = Float64.(collect(lons)),
            ys = Float64.(collect(lats)),
            fields = [Float64.(F) for F in fields],
            xlabel = "Longitude (°E)",
            ylabel = "Latitude (°N)",
            shelf_frame = nothing
        )
    end
end

"""
    _overlay_coastline!(ax, coastline; fill_color, stroke_color, linewidth, shelf_frame)

Draw filled land polygons for closed coastline rings and crisp shoreline strokes for
open continental arcs on Makie axis `ax`. When `shelf_frame` is provided as
`(lon0, lat0, rot_deg)`, coordinates are projected into the shelf-aligned km frame.
"""
function _overlay_coastline!(
    ax,
    coastline;
    fill_color = :white,
    stroke_color = :gray20,
    linewidth::Real = 1.0,
    shelf_frame::Union{Nothing, Tuple{Float64, Float64, Float64}} = nothing
)
    isnothing(coastline) && return nothing
    for p in coastline
        length(p.lons) >= 2 || continue
        pts = if isnothing(shelf_frame)
            Point2f.(Float64.(p.lons), Float64.(p.lats))
        else
            lon0, lat0, rot_deg = shelf_frame
            [let (xk, yk) = _lonlat_to_shelf_km(p.lons[k], p.lats[k], lon0, lat0, rot_deg)
                 Point2f(xk, yk)
             end for k in eachindex(p.lons)]
        end
        is_closed = length(p.lons) >= 3 &&
                    p.lons[begin] == p.lons[end] &&
                    p.lats[begin] == p.lats[end]
        if is_closed
            poly!(ax, pts; color = fill_color, strokecolor = stroke_color, strokewidth = linewidth)
        else
            lines!(ax, pts; color = stroke_color, linewidth = linewidth)
        end
    end
    return nothing
end

"""
    _focus_wet_limits!(ax, lons, lats, field; margin = 0.3)

Set `xlims!` and `ylims!` tightly around the finite (wet) cells of `field` with a small
degree margin so maps focus on the ocean area of interest rather than empty margins.
"""
function _focus_wet_limits!(ax, lons, lats, field; margin::Real = 0.3)
    wet_ij = findall(isfinite, field)
    if !isempty(wet_ij)
        wet_lons = [lons[idx[1]] for idx in wet_ij]
        wet_lats = [lats[idx[2]] for idx in wet_ij]
        xlims!(ax, minimum(wet_lons) - margin, maximum(wet_lons) + margin)
        ylims!(ax, minimum(wet_lats) - margin, maximum(wet_lats) + margin)
    else
        xlims!(ax, minimum(lons), maximum(lons))
        ylims!(ax, minimum(lats), maximum(lats))
    end
    return nothing
end

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
    coastline = nothing,
    rotation_deg::Real = 0.0,
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

    use_rot = abs(Float64(rotation_deg)) > 1e-3
    if use_rot
        res = _resample_to_shelf_grid(lons, lats, spd_k, spd_k, u_k, v_k; rot_deg = rotation_deg)
        xs, ys = res.xs, res.ys
        spd_plot = res.fields[1]
        c_r, s_r = cosd(res.rot_deg), sind(res.rot_deg)
        u_plot = res.fields[2] .* c_r .+ res.fields[3] .* s_r
        v_plot = -res.fields[2] .* s_r .+ res.fields[3] .* c_r
        xlabel_txt = "Along-Shelf Distance (km, $(round(Int, res.rot_deg))° CCW from E)"
        ylabel_txt = "Cross-Shelf Distance (km)"
        shelf_frame = (res.lon0, res.lat0, res.rot_deg)
    else
        xs, ys = Float64.(collect(lons)), Float64.(collect(lats))
        spd_plot, u_plot, v_plot = spd_k, u_k, v_k
        xlabel_txt = "Longitude (°E)"
        ylabel_txt = "Latitude (°N)"
        shelf_frame = nothing
    end

    fig = Figure(size = (1020, 640), fontsize = 13)
    ax = Axis(
        fig[1, 1],
        title = "$title\n[Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = xlabel_txt,
        ylabel = ylabel_txt,
        aspect = DataAspect()
    )

    finite_spd = sort(filter(isfinite, vec(spd_plot)))
    max_spd = isempty(finite_spd) ? 1.0 : max(1.0, finite_spd[end])
    p90_spd = isempty(finite_spd) ? 1.0 : max(1.0, finite_spd[clamp(round(Int, 0.90 * length(finite_spd)), 1, length(finite_spd))])
    hm = heatmap!(ax, xs, ys, spd_plot, colormap = colormap;
                  colorrange = (0, max_spd), nan_color = :transparent)
    Colorbar(fig[1, 2], hm, label = "Current Speed (cm/s)")

    # Quiver overlay: ~26 vectors along the longer axis, scaled by p90 speed so shelf
    # currents remain clearly visible without outliers shrinking the vector field.
    stride = max(1, quiver_stride, round(Int, max(length(xs), length(ys)) / 26))
    ix = 1:stride:length(xs)
    iy = 1:stride:length(ys)
    dx_q = length(ix) > 1 ? abs(xs[ix[end]] - xs[ix[1]]) / (length(ix) - 1) : 1.0
    dy_q = length(iy) > 1 ? abs(ys[iy[end]] - ys[iy[1]]) / (length(iy) - 1) : 1.0
    arrow_scale = 1.45 * max(dx_q, dy_q) / p90_spd
    pts = Point2f[]
    dirs = Vec2f[]
    for i in ix, j in iy
        uij = u_plot[i, j]
        vij = v_plot[i, j]
        if isfinite(uij) && isfinite(vij)
            mag = hypot(uij, vij)
            cap = mag > 2.2 * p90_spd ? (2.2 * p90_spd / mag) : 1.0
            push!(pts, Point2f(xs[i], ys[j]))
            push!(dirs, Vec2f(uij * cap * arrow_scale, vij * cap * arrow_scale))
        end
    end
    if !isempty(pts)
        arrows2d!(
            ax,
            pts,
            dirs;
            tipwidth = 9,
            tiplength = 11,
            shaftwidth = 1.8,
            color = (:white, 0.92)
        )
    end

    # Bathymetry contours (geographic mode)
    if !use_rot && !isnothing(bathymetry_data) && hasproperty(bathymetry_data, :elevation)
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

    _overlay_coastline!(ax, coastline; shelf_frame = shelf_frame)
    if use_rot
        xlims!(ax, xs[1], xs[end])
        ylims!(ax, ys[1], ys[end])
    else
        _focus_wet_limits!(ax, xs, ys, spd_plot)
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
    coastline = nothing,
    rotation_deg::Real = 0.0,
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
    mp = _prepare_map_fields(lons, lats, T_k, T_k, S_k; rotation_deg = rotation_deg)
    T_p, S_p = mp.fields[1], mp.fields[2]

    fig = Figure(size = (1200, 480), fontsize = 12)

    ax1 = Axis(fig[1, 1],
        title = "Surface Temperature (°C)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm1 = heatmap!(ax1, mp.xs, mp.ys, T_p, colormap = colormap_temp,
                   colorrange = _finite_colorrange(T_p; fallback = (0.0, 15.0)),
                   nan_color = :transparent)
    _overlay_coastline!(ax1, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax1, mp.xs, mp.ys, T_p)
    else
        xlims!(ax1, mp.xs[1], mp.xs[end]); ylims!(ax1, mp.ys[1], mp.ys[end])
    end
    Colorbar(fig[1, 2], hm1, label = "°C")

    ax2 = Axis(fig[1, 3],
        title = "Surface Salinity (PSU)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm2 = heatmap!(ax2, mp.xs, mp.ys, S_p, colormap = colormap_salt,
                   colorrange = _finite_colorrange(S_p; fallback = (30.0, 35.0)),
                   nan_color = :transparent)
    _overlay_coastline!(ax2, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax2, mp.xs, mp.ys, T_p)
    else
        xlims!(ax2, mp.xs[1], mp.xs[end]); ylims!(ax2, mp.ys[1], mp.ys[end])
    end
    Colorbar(fig[1, 4], hm2, label = "PSU")

    fig[0, :] = Label(fig, "Hydrodynamic Tracers $(isempty(title) ? "" : "- $title")\n[Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
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
    coastline = nothing,
    rotation_deg::Real = 0.0,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_stratification.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    N2_k = Float64.(hydro.stratification[:, :, k] .* 1e4)  # convert to 10^-4 s^-2
    mp = _prepare_map_fields(lons, lats, N2_k, N2_k; rotation_deg = rotation_deg)
    N2_p = mp.fields[1]

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Buoyancy Frequency N² — [Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm = heatmap!(ax, mp.xs, mp.ys, N2_p, colormap = :ice;
                   colorrange = _finite_colorrange(N2_p; fallback = (0.0, 5.0)),
                   nan_color = :transparent)
    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax, mp.xs, mp.ys, N2_p)
    else
        xlims!(ax, mp.xs[1], mp.xs[end]); ylims!(ax, mp.ys[1], mp.ys[end])
    end
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
    coastline = nothing,
    rotation_deg::Real = 0.0,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_diffusion.png"
)
    lons = hydro.lons
    lats = hydro.lats
    depths = hydro.depths
    k = _surface_depth_index(depths, depth_level)
    depth_m = depths[k]
    t_hr = hydro.time_seconds / 3600.0

    κ_k = Float64.(hydro.diffusion[:, :, k])
    κ_log = log10.(clamp.(κ_k, 1e-8, 1e-2))
    mp = _prepare_map_fields(lons, lats, hydro.u[:, :, k], κ_log; rotation_deg = rotation_deg)
    κ_p = mp.fields[1]

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Vertical Eddy Diffusivity κ_v — [Depth: $(round(depth_m, digits=1)) m | Time: $(round(t_hr, digits=1)) h]",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm = heatmap!(ax, mp.xs, mp.ys, κ_p, colormap = :turbid;
                   colorrange = _finite_colorrange(κ_p; fallback = (-8.0, -2.0)),
                   nan_color = :transparent)
    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax, lons, lats, hydro.u[:, :, k])
    else
        xlims!(ax, mp.xs[1], mp.xs[end]); ylims!(ax, mp.ys[1], mp.ys[end])
    end
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
    section_type::Symbol = :cross_shelf,
    transect_endpoints::Union{Nothing, Tuple{Tuple{Real, Real}, Tuple{Real, Real}}} = nothing,
    n_samples::Int = 160,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_section.png"
)
    lons = Float64.(collect(hydro.lons))
    lats = Float64.(collect(hydro.lats))
    depths = Float64.(collect(hydro.depths))
    nz = length(depths)

    field = variable === :temperature ? hydro.temperature : hydro.salinity
    var_label = variable === :temperature ? "Temperature" : "Salinity"
    bathy = hydro.bathymetry

    if section_type === :cross_shelf || section_type === :along_shelf || !isnothing(transect_endpoints)
        # Default Scotian Shelf transects:
        #   :cross_shelf -> NW-to-SE perpendicular to the Nova Scotia coastline (Halifax -> shelf break)
        #   :along_shelf -> SW-to-NE parallel to the shelf axis
        (lon_a, lat_a), (lon_b, lat_b) = if !isnothing(transect_endpoints)
            transect_endpoints
        elseif section_type === :along_shelf
            ((-65.8, 42.7), (-58.8, 45.7))
        else
            ((-63.6, 44.6), (-60.8, 42.2))
        end
        t_lons = collect(range(Float64(lon_a), Float64(lon_b); length = n_samples))
        t_lats = collect(range(Float64(lat_a), Float64(lat_b); length = n_samples))
        mid_lat = 0.5 * (lat_a + lat_b)
        total_km = hypot((lon_b - lon_a) * 111.320 * cosd(mid_lat), (lat_b - lat_a) * 111.132)
        x_axis = collect(range(0.0, total_km; length = n_samples))
        sec = fill(NaN, n_samples, nz)
        for iz in 1:nz
            slice_z = Float64.(field[:, :, iz])
            for ix in 1:n_samples
                val = _bilinear_2d(lons, lats, slice_z, t_lons[ix], t_lats[ix])
                b_val = _bilinear_2d(lons, lats, bathy, t_lons[ix], t_lats[ix])
                if isfinite(val) && (!isfinite(b_val) || depths[iz] >= b_val)
                    sec[ix, iz] = val
                end
            end
        end
        x_label = section_type === :along_shelf ?
                  "Along-Shelf Distance (km, SW → NE)" :
                  "Cross-Shelf Distance Offshore (km, Coast → Slope)"
        title = "$(var_label) $(section_type === :along_shelf ? "Along-Shelf" : "Cross-Shelf") Section"
    elseif section_type === :lat
        j = argmin(abs.(lats .- coordinate))
        sec = Float64.(field[:, j, :])
        x_axis = lons
        x_label = "Longitude (°E)"
        title = "$(var_label) Cross-Section at Latitude $(round(coordinate, digits=1))°N"
        for ix in 1:length(x_axis), iz in 1:nz
            b_val = bathy[ix, j]
            if isfinite(b_val) && depths[iz] < b_val
                sec[ix, iz] = NaN
            end
        end
    else
        i = argmin(abs.(lons .- coordinate))
        sec = Float64.(field[i, :, :])
        x_axis = lats
        x_label = "Latitude (°N)"
        title = "$(var_label) Cross-Section at Longitude $(round(coordinate, digits=1))°E"
        for ix in 1:length(x_axis), iz in 1:nz
            b_val = bathy[i, ix]
            if isfinite(b_val) && depths[iz] < b_val
                sec[ix, iz] = NaN
            end
        end
    end

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = title,
        xlabel = x_label,
        ylabel = "Depth (m)"
    )
    hm = heatmap!(ax, x_axis, depths, sec;
                  colormap = variable === :temperature ? :thermal : :haline,
                  colorrange = _finite_colorrange(sec; fallback = (0.0, 15.0)),
                  nan_color = :transparent)
    Colorbar(fig[1, 2], hm, label = variable === :temperature ? "°C" : "PSU")

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

function _nearest_wet_station(lons, lats, field_surf, station)
    target_lon = isnothing(station) ? 0.5 * (lons[1] + lons[end]) : Float64(station[1])
    target_lat = isnothing(station) ? 0.5 * (lats[1] + lats[end]) : Float64(station[2])
    i_st = argmin(abs.(lons .- target_lon))
    j_st = argmin(abs.(lats .- target_lat))
    if isfinite(field_surf[i_st, j_st])
        return (i_st, j_st)
    end
    best_i, best_j, best_d2 = i_st, j_st, Inf
    for j in eachindex(lats), i in eachindex(lons)
        isfinite(field_surf[i, j]) || continue
        d2 = (lons[i] - target_lon)^2 + (lats[j] - target_lat)^2
        if d2 < best_d2
            best_d2 = d2
            best_i, best_j = i, j
        end
    end
    return (best_i, best_j)
end

"""
    plot_hydrodynamic_timeseries(
        hydro::NamedTuple;
        store = nothing,
        station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        variable::Symbol = :temperature,
        max_samples::Int = 200,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
    ) -> Figure

Plot surface and near-bed Eulerian time series at a representative wet shelf station.
When `store` is supplied, reads the station column across the full simulation record.
"""
function plot_hydrodynamic_timeseries(
    hydro::NamedTuple;
    store = nothing,
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    variable::Symbol = :temperature,
    max_samples::Int = 200,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_timeseries.png"
)
    lons = Float64.(collect(hydro.lons))
    lats = Float64.(collect(hydro.lats))
    depths = Float64.(collect(hydro.depths))
    k_surf = argmax(depths)

    field = variable === :temperature ? hydro.temperature : hydro.salinity
    i_st, j_st = _nearest_wet_station(lons, lats, field[:, :, k_surf], station)

    # Find deepest wet level at this station
    k_bed = k_surf
    for k in sortperm(depths)
        if isfinite(field[i_st, j_st, k])
            k_bed = k
            break
        end
    end

    var_name = variable === :temperature ? "T" : "S"
    times_days = Float64[]
    ts_surf = Float64[]
    ts_bed = Float64[]
    if !isnothing(store) && (var_name in store.variables)
        nt = length(store.times)
        keep = nt <= max_samples ? collect(1:nt) :
               sort(unique(round.(Int, range(1, nt; length = max_samples))))
        t0 = store.times[1]
        arr = store.dataset.variables[var_name].data
        for m in keep
            push!(times_days, (Float64(store.times[m]) - t0) / 86400.0)
            push!(ts_surf, Float64(arr[i_st, j_st, k_surf, m]))
            push!(ts_bed, Float64(arr[i_st, j_st, k_bed, m]))
        end
    else
        push!(times_days, hydro.time_seconds / 86400.0)
        push!(ts_surf, Float64(field[i_st, j_st, k_surf]))
        push!(ts_bed, Float64(field[i_st, j_st, k_bed]))
    end

    fig = Figure(size = (760, 380), fontsize = 12)
    ax = Axis(fig[1, 1],
        title = "Eulerian $(variable) Time Series at Shelf Station ($(round(lons[i_st], digits=2))°E, $(round(lats[j_st], digits=2))°N)",
        xlabel = "Elapsed Time (days)",
        ylabel = variable === :temperature ? "Temperature (°C)" : "Salinity (PSU)"
    )
    lines!(ax, times_days, ts_surf; color = :crimson, linewidth = 2,
           label = "Surface ($(round(depths[k_surf], digits=1)) m)")
    if k_bed != k_surf
        lines!(ax, times_days, ts_bed; color = :steelblue, linewidth = 2,
               label = "Near-bed ($(round(depths[k_bed], digits=1)) m)")
    end
    if length(times_days) == 1
        scatter!(ax, times_days, ts_surf; color = :crimson, markersize = 6)
    end
    axislegend(ax; position = :rt)

    if !isnothing(output_path)
        mkpath(dirname(output_path))
        save(output_path, fig)
    end

    return fig
end

"""
    plot_hydrodynamic_hovmoller(
        hydro::NamedTuple;
        store = nothing,
        station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
        variable::Symbol = :temperature,
        max_samples::Int = 200,
        output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_depth_time_profile.png"
    ) -> Figure

Plot a vertical depth–time profile (water-column evolution over time) at a wet shelf station.
"""
function plot_hydrodynamic_hovmoller(
    hydro::NamedTuple;
    store = nothing,
    station::Union{Nothing, Tuple{<:Real, <:Real}} = nothing,
    variable::Symbol = :temperature,
    max_samples::Int = 200,
    output_path::Union{Nothing, AbstractString} = "outputs/hydrodynamic_depth_time_profile.png"
)
    lons = Float64.(collect(hydro.lons))
    lats = Float64.(collect(hydro.lats))
    depths = Float64.(collect(hydro.depths))
    k_surf = argmax(depths)
    nz = length(depths)

    field = variable === :temperature ? hydro.temperature : hydro.salinity
    i_st, j_st = _nearest_wet_station(lons, lats, field[:, :, k_surf], station)

    var_name = variable === :temperature ? "T" : "S"
    if !isnothing(store) && (var_name in store.variables)
        nt = length(store.times)
        keep = nt <= max_samples ? collect(1:nt) :
               sort(unique(round.(Int, range(1, nt; length = max_samples))))
        t0 = store.times[1]
        times_days = [(Float64(store.times[m]) - t0) / 86400.0 for m in keep]
        arr = store.dataset.variables[var_name].data
        sec = Matrix{Float64}(undef, length(keep), nz)
        for (idx, m) in enumerate(keep)
            sec[idx, :] .= Float64.(arr[i_st, j_st, :, m])
        end
    else
        times_days = [hydro.time_seconds / 86400.0]
        sec = reshape(Float64.(field[i_st, j_st, :]), 1, nz)
    end

    # Restrict y-axis to wet water-column depths at this station
    wet_k = [k for k in 1:nz if any(isfinite, sec[:, k])]
    min_z = isempty(wet_k) ? minimum(depths) : minimum(depths[wet_k])

    fig = Figure(size = (820, 440), fontsize = 12)
    ax = Axis(fig[1, 1],
        title = "$(titlecase(string(variable))) Depth–Time Profile at Station ($(round(lons[i_st], digits=2))°E, $(round(lats[j_st], digits=2))°N)",
        xlabel = "Elapsed Time (days)",
        ylabel = "Depth (m)"
    )
    hm = heatmap!(ax, times_days, depths, sec;
                  colormap = variable === :temperature ? :thermal : :haline,
                  colorrange = _finite_colorrange(sec; fallback = (0.0, 15.0)),
                  nan_color = :transparent)
    ylims!(ax, min_z, 0.0)
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
    coastline = nothing,
    rotation_deg::Real = 0.0,
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

    mp = _prepare_map_fields(lons, lats, T_bottom, T_bottom; rotation_deg = rotation_deg)
    T_p = mp.fields[1]

    fig = Figure(size = (980, 480), fontsize = 13)
    ax = Axis(fig[1, 1],
        title = "Bottom Temperature (°C)",
        xlabel = mp.xlabel,
        ylabel = mp.ylabel,
        aspect = DataAspect()
    )
    hm = heatmap!(ax, mp.xs, mp.ys, T_p, colormap = :thermal;
                   colorrange = _finite_colorrange(T_p; fallback = (0.0, 15.0)),
                   nan_color = :transparent)
    _overlay_coastline!(ax, coastline; shelf_frame = mp.shelf_frame)
    if isnothing(mp.shelf_frame)
        _focus_wet_limits!(ax, mp.xs, mp.ys, T_p)
    else
        xlims!(ax, mp.xs[1], mp.xs[end]); ylims!(ax, mp.ys[1], mp.ys[end])
    end
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



