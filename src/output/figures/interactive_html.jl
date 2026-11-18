"""
    interactive_html.jl

Interactive Leaflet dashboard engine.

The standalone HTML dashboard and the JavaScript it embeds. One file, because
the embedded script calls the helpers and the helpers build the script;
separating them would split a template from the code that fills it.
"""

using CairoMakie
import CairoMakie: record
using Zarr

# Format a real scalar for RFC 8259 JSON, mapping non-finite values to `null`.
_json_num(x::Real; digits::Int = 4) = isfinite(x) ? string(round(x; digits = digits)) : "null"

"""
    export_interactive_tracks_html(
        output_path::AbstractString;
        trajectories::Union{Nothing, NamedTuple} = nothing,
        runs_data::Union{Nothing, AbstractDict} = nothing,
        hydrodynamics::Union{Nothing, NamedTuple, AbstractDict, Any} = nothing,
        gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
        connectivity::Union{Nothing, NamedTuple} = nothing,
        recruitment_metrics::Union{Nothing, NamedTuple} = nothing,
        bathymetry::Union{Nothing, NamedTuple} = nothing,
        strata_definitions::Union{Nothing, AbstractVector} = nothing,
        title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map",
        max_tracks::Int = 0,
        max_points_per_track::Int = 0
    ) -> String

Generate an interactive HTML5 + Leaflet.js dashboard paired with a
zlib-compressed binary sidecar (`<stem>.data.zz`) containing the 4D Lagrangian
particle tracks, developmental stages, DVM behaviors, settlement density, and
regional Crab Fishing Area (CFA) connectivity. Eulerian hydrodynamic grids are
intentionally omitted from the dashboard payload (handled by TethysFlow).

When served over HTTP/HTTPS (including GitHub Pages or a local static server),
the HTML fetches `<stem>.data.zz` automatically and decompresses it in the
browser via `DecompressionStream("deflate")`. When opened directly over
`file://` where browser CORS blocks `fetch()`, a file-picker overlay allows
selecting `<stem>.data.zz` locally.

# Mathematical & Behavioral Visualization
- **4D Trajectories**: \$\\mathbf{x}_p(t) = (\\\\\\\\lambda_p(t), \\phi_p(t), z_p(t))\$ with
  temperature exposure \$T_p(t)\$ and cumulative degree-days \$DD_p(t)\$.
- **Developmental Molting Coloration**: Zoea I \$\\to\$ IV, Megalopa, Instar I (settled).
- **Playback Animation**: Interactive time-stepping simulation scrubber and velocity trails.

# Inputs
- `output_path::AbstractString`: Destination path (e.g. `"outputs/interactive_larval_tracks.html"`).
  A companion sidecar `<stem>.data.zz` is written alongside it.
- `trajectories::NamedTuple`: Particle tracking trajectory NamedTuple.
- `runs_data`: Optional dictionary mapping run names to data bundles.
- `hydrodynamics`: Retained as a no-op keyword for signature compatibility; ignored.
- `gridded_dispersal`: Optional 2D spatial dispersal fields.
- `connectivity`: Optional connectivity transition matrix.
- `recruitment_metrics`: Optional summary recruitment metrics.
- `bathymetry`: Retained as a no-op keyword for signature compatibility; ignored.
- `strata_definitions`: Optional spatial management boundaries.
- `title::AbstractString`: Dashboard header title.
- `max_tracks::Int`: Optional upper bound on exported tracks (`<= 0` exports all tracks).
- `max_points_per_track::Int`: Optional upper bound on time points per track
  (`<= 0` exports all time points).

# Outputs
- `String`: Path to the generated interactive HTML file.
"""
function export_interactive_tracks_html(
    output_path::AbstractString;
    trajectories::Union{Nothing, NamedTuple} = nothing,
    runs_data::Union{Nothing, AbstractDict} = nothing,
    hydrodynamics::Union{Nothing, NamedTuple, AbstractDict, Any} = nothing,
    gridded_dispersal::Union{Nothing, NamedTuple} = nothing,
    connectivity::Union{Nothing, NamedTuple} = nothing,
    recruitment_metrics::Union{Nothing, NamedTuple} = nothing,
    bathymetry::Union{Nothing, NamedTuple} = nothing,
    strata_definitions::Union{Nothing, AbstractVector} = nothing,
    title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map",
    max_tracks::Int = 0,
    max_points_per_track::Int = 0
)::String
    out_dir = dirname(output_path)
    if !isempty(out_dir) && !isdir(out_dir)
        mkpath(out_dir)
    end

    # 1. Structure runs dictionary
    runs_dict = Dict{String, Any}()
    if !isnothing(runs_data) && !isempty(runs_data)
        for (k, v) in runs_data
            key_str = string(k)
            t_obj = hasproperty(v, :trajectories) ? v.trajectories : (v isa NamedTuple && hasproperty(v, :lons) ? v : nothing)
            if !isnothing(t_obj)
                runs_dict[key_str] = (
                    trajectories = t_obj,
                    gridded_dispersal = hasproperty(v, :gridded_dispersal) ? v.gridded_dispersal : (hasproperty(v, :empirical_movement) ? v.empirical_movement : gridded_dispersal),
                    connectivity = hasproperty(v, :connectivity) ? v.connectivity : connectivity,
                    recruitment_metrics = hasproperty(v, :recruitment_metrics) ? v.recruitment_metrics : recruitment_metrics,
                    hydrodynamics = hasproperty(v, :hydrodynamics) ? v.hydrodynamics : (hasproperty(v, :model) ? v.model : hydrodynamics)
                )
            end
        end
    elseif !isnothing(trajectories)
        runs_dict["Active Simulation"] = (
            trajectories = trajectories,
            gridded_dispersal = gridded_dispersal,
            connectivity = connectivity,
            recruitment_metrics = recruitment_metrics,
            hydrodynamics = hydrodynamics
        )
    else
        error("Either `trajectories` or `runs_data` must be provided to `export_interactive_tracks_html`.")
    end

    # 2. Extract Management Strata and Centroids
    strata_list = if isnothing(strata_definitions)
        auto_polys = load_cfa_polygons("inputs")
        if !isempty(auto_polys)
            [
                (
                    name = string(s.name),
                    color = s.color,
                    is_polygon = true,
                    polygon = s.coordinates,
                    lon_min = s.lon[1],
                    lon_max = s.lon[2],
                    lat_min = s.lat[1],
                    lat_max = s.lat[2],
                    centroid_lat = mean([pt[1] for pt in s.coordinates]),
                    centroid_lon = mean([pt[2] for pt in s.coordinates])
                ) for s in auto_polys
            ]
        else
            [
                (name = "CFA 20-22 (Eastern NS)", color = "#3B82F6", is_polygon = false, polygon = Any[], lon_min = -62.0, lon_max = -57.0, lat_min = 44.5, lat_max = 47.5, centroid_lat = 46.0, centroid_lon = -59.5),
                (name = "CFA 23-24 (Middle Shelf)", color = "#10B981", is_polygon = false, polygon = Any[], lon_min = -64.5, lon_max = -60.0, lat_min = 43.0, lat_max = 45.5, centroid_lat = 44.25, centroid_lon = -62.25),
                (name = "CFA 4X (Southwest NS)", color = "#F59E0B", is_polygon = false, polygon = Any[], lon_min = -67.5, lon_max = -63.5, lat_min = 42.0, lat_max = 44.5, centroid_lat = 43.25, centroid_lon = -65.5),
                (name = "Offshore / Slope", color = "#8B5CF6", is_polygon = false, polygon = Any[], lon_min = -68.0, lon_max = -57.0, lat_min = 41.5, lat_max = 43.0, centroid_lat = 42.25, centroid_lon = -62.5)
            ]
        end
    else
        [
            (
                name = string(s.name),
                color = hasproperty(s, :color) ? string(s.color) : "#3B82F6",
                is_polygon = hasproperty(s, :coordinates) && length(s.coordinates) >= 3,
                polygon = hasproperty(s, :coordinates) ? s.coordinates :
                          (hasproperty(s, :lons) && hasproperty(s, :lats) ? [[s.lats[i], s.lons[i]] for i in 1:length(s.lons)] : Any[]),
                lon_min = hasproperty(s, :lon) ? Float64(s.lon[1]) : (hasproperty(s, :lon_min) ? Float64(s.lon_min) : -68.0),
                lon_max = hasproperty(s, :lon) ? Float64(s.lon[2]) : (hasproperty(s, :lon_max) ? Float64(s.lon_max) : -57.0),
                lat_min = hasproperty(s, :lat) ? Float64(s.lat[1]) : (hasproperty(s, :lat_min) ? Float64(s.lat_min) : 42.0),
                lat_max = hasproperty(s, :lat) ? Float64(s.lat[2]) : (hasproperty(s, :lat_max) ? Float64(s.lat_max) : 47.0),
                centroid_lat = hasproperty(s, :coordinates) && !isempty(s.coordinates) ? mean([pt[1] for pt in s.coordinates]) : (hasproperty(s, :lat) ? (s.lat[1]+s.lat[2])/2.0 : 44.0),
                centroid_lon = hasproperty(s, :coordinates) && !isempty(s.coordinates) ? mean([pt[2] for pt in s.coordinates]) : (hasproperty(s, :lon) ? (s.lon[1]+s.lon[2])/2.0 : -63.0)
            ) for s in strata_definitions
        ]
    end

    strata_json_entries = String[]
    for s in strata_list
        poly_str = if s.is_polygon && !isempty(s.polygon)
            "[" * join(["[$(pt[1]),$(pt[2])]" for pt in s.polygon], ",") * "]"
        else
            "[]"
        end
        s_entry = """{
            "name": "$(s.name)",
            "color": "$(s.color)",
            "is_polygon": $(s.is_polygon),
            "polygon": $(poly_str),
            "lon_min": $(s.lon_min),
            "lon_max": $(s.lon_max),
            "lat_min": $(s.lat_min),
            "lat_max": $(s.lat_max),
            "centroid_lat": $(s.centroid_lat),
            "centroid_lon": $(s.centroid_lon)
        }"""
        push!(strata_json_entries, s_entry)
    end
    strata_json = "[" * join(strata_json_entries, ",\n") * "]"

    # 3. Serialize each run
    run_json_blocks = String[]
    for (run_key, run_data) in runs_dict
        trajs = canonicalize_trajectories(run_data.trajectories)
        lons = trajs.lons
        lats = trajs.lats
        depths = trajs.depths
        stages = trajs.stages
        times = trajs.times
        n_p, n_t = size(lons)

        temps = trajs.temperatures
        dds = trajs.degree_days_timeseries
        survs = trajs.survival_probability

        statuses = [string(s) for s in trajs.settlement_status]
        alive = trajs.alive
        ages = trajs.settlement_age
        ids = trajs.ids

        # Optional uniform decimation when max_tracks > 0 or max_points_per_track > 0.
        # Defaults (<= 0) retain all tracks and time steps in the compressed sidecar;
        # summary statistics are always computed from full-resolution arrays.
        track_idx = if max_tracks > 0 && n_p > max_tracks
            unique(round.(Int, range(1, n_p, length = max_tracks)))
        else
            collect(1:n_p)
        end
        time_idx = if max_points_per_track > 0 && n_t > max_points_per_track
            unique(round.(Int, range(1, n_t, length = max_points_per_track)))
        else
            collect(1:n_t)
        end
        n_tracks_shown = length(track_idx)
        n_points_shown = length(time_idx)
        t_days = [round(times[t] / 86400.0, digits = 3) for t in time_idx]
        decimated = (n_tracks_shown < n_p) || (n_points_shown < n_t)
        if decimated
            @info "export_interactive_tracks_html: decimated payload to " *
                  "$n_tracks_shown of $n_p tracks and $n_points_shown of $n_t time points. " *
                  "Statistics are computed from full-resolution arrays and are unaffected."
        end

        p_entries = String[]
        for p in track_idx
            p_lons = [_json_num(lons[p, t]; digits = 4) for t in time_idx]
            p_lats = [_json_num(lats[p, t]; digits = 4) for t in time_idx]
            p_depths = [_json_num(depths[p, t]; digits = 0) for t in time_idx]
            p_temps = [_json_num(temps[p, t]; digits = 1) for t in time_idx]
            p_dds = [_json_num(dds[p, t]; digits = 1) for t in time_idx]
            p_survs = [_json_num(survs[p, t]; digits = 3) for t in time_idx]
            p_stages = [string(stages[p, t]) for t in time_idx]

            push!(p_entries, """{
                "id": $(ids[p]),
                "status": "$(statuses[p])",
                "alive": $(alive[p]),
                "settlement_age_days": $(_json_num(ages[p] / 86400.0; digits = 2)),
                "lons": [$(join(p_lons, ","))],
                "lats": [$(join(p_lats, ","))],
                "depths": [$(join(p_depths, ","))],
                "temps": [$(join(p_temps, ","))],
                "dds": [$(join(p_dds, ","))],
                "survs": [$(join(p_survs, ","))],
                "stages": [$(join(["\"$s\"" for s in p_stages], ","))]
            }""")
        end
        particles_json = "[" * join(p_entries, ",\n") * "]"

        # Statistics
        n_settled = count(==("settled_successful"), statuses)
        n_alive = count(identity, alive)
        settle_pct = n_p > 0 ? round(100.0 * n_settled / n_p, digits = 1) : 0.0
        survival_pct = n_p > 0 ? round(100.0 * n_alive / n_p, digits = 1) : 0.0
        finite_ages = filter(isfinite, [Float64(a) / 86400.0 for a in ages])
        mean_pld = isempty(finite_ages) ? 0.0 : round(mean(finite_ages), digits = 1)
        finite_temps = filter(isfinite, vec(temps))
        mean_temp_val = isempty(finite_temps) ? 0.0 : round(mean(finite_temps), digits = 2)

        # Gridded settlement density, velocity quivers, and thermal field
        g_disp = run_data.gridded_dispersal
        dens_entries = String[]
        quiver_entries = String[]
        therm_entries = String[]

        if !isnothing(g_disp) && hasproperty(g_disp, :lon_centers) && hasproperty(g_disp, :lat_centers)
            lons_c = g_disp.lon_centers
            lats_c = g_disp.lat_centers
            nx_c = length(lons_c)
            ny_c = length(lats_c)
            step_stride = max(1, round(Int, nx_c / 25))

            has_dens = hasproperty(g_disp, :density) || hasproperty(g_disp, :settlement_density)
            dens_mat = hasproperty(g_disp, :density) ? g_disp.density : (hasproperty(g_disp, :settlement_density) ? g_disp.settlement_density : nothing)

            has_u = hasproperty(g_disp, :u_mean)
            has_v = hasproperty(g_disp, :v_mean)
            has_diff = hasproperty(g_disp, :diffusivity)
            has_t = hasproperty(g_disp, :mean_exposure_temperature)
            has_dd = hasproperty(g_disp, :mean_degree_days)

            for i in 1:step_stride:nx_c, j in 1:step_stride:ny_c
                lon_val = round(lons_c[i], digits = 4)
                lat_val = round(lats_c[j], digits = 4)

                if !isnothing(dens_mat) && !isnan(dens_mat[i, j]) && dens_mat[i, j] > 1e-4
                    d_val = round(dens_mat[i, j], digits = 5)
                    push!(dens_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"density\":$(d_val)}")
                end

                if has_u && has_v && !isnan(g_disp.u_mean[i, j]) && !isnan(g_disp.v_mean[i, j])
                    u_cm = round(g_disp.u_mean[i, j] * 100.0, digits = 2)
                    v_cm = round(g_disp.v_mean[i, j] * 100.0, digits = 2)
                    spd = round(hypot(u_cm, v_cm), digits = 2)
                    diff_val = has_diff && !isnan(g_disp.diffusivity[i, j]) ? round(g_disp.diffusivity[i, j], digits = 1) : 10.0
                    if spd > 0.1
                        push!(quiver_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"u\":$(u_cm),\"v\":$(v_cm),\"speed\":$(spd),\"diff\":$(diff_val)}")
                    end
                end

                if has_t && !isnan(g_disp.mean_exposure_temperature[i, j])
                    t_val = round(g_disp.mean_exposure_temperature[i, j], digits = 2)
                    dd_val = has_dd && !isnan(g_disp.mean_degree_days[i, j]) ? round(g_disp.mean_degree_days[i, j], digits = 1) : 40.0
                    push!(therm_entries, "{\"lon\":$(lon_val),\"lat\":$(lat_val),\"temp\":$(t_val),\"dd\":$(dd_val)}")
                end
            end
        end

        # Demographic connectivity flows
        conn = run_data.connectivity
        flow_entries = String[]
        if !isnothing(conn) && hasproperty(conn, :matrix) && hasproperty(conn, :strata_names)
            mat = conn.matrix
            s_names = conn.strata_names
            n_st = length(s_names)

            c_map = Dict(s.name => (s.centroid_lat, s.centroid_lon) for s in strata_list)

            for i in 1:n_st, j in 1:n_st
                p_ij = Float64(mat[i, j])
                if p_ij > 0.005
                    src_n = s_names[i]
                    dst_n = s_names[j]
                    src_c = get(c_map, src_n, (44.0, -63.0))
                    dst_c = get(c_map, dst_n, (44.0, -63.0))
                    push!(flow_entries, """{
                        "src": "$(src_n)",
                        "dst": "$(dst_n)",
                        "src_lat": $(round(src_c[1], digits=4)),
                        "src_lon": $(round(src_c[2], digits=4)),
                        "dst_lat": $(round(dst_c[1], digits=4)),
                        "dst_lon": $(round(dst_c[2], digits=4)),
                        "prob": $(round(p_ij, digits=3))
                    }""")
                end
            end
        end

        run_json = """
        "$(run_key)": {
            "name": "$(run_key)",
            "times": [$(join(t_days, ","))],
            "stats": {
                "n_particles": $(n_p),
                "n_settled": $(n_settled),
                "n_alive": $(n_alive),
                "settle_pct": $(settle_pct),
                "survival_pct": $(survival_pct),
                "mean_pld_days": $(mean_pld),
                "mean_temp": $(mean_temp_val)
            },
            "particles": $(particles_json),
            "settlement_density": [$(join(dens_entries, ","))],
            "velocity_quivers": [$(join(quiver_entries, ","))],
            "thermal_field": [$(join(therm_entries, ","))],
            "connectivity_flows": [$(join(flow_entries, ","))]
        }
        """
        push!(run_json_blocks, run_json)
    end

    runs_master_json = "{\n" * join(run_json_blocks, ",\n") * "\n}"
    default_run_key = first(keys(runs_dict))
    sidecar_path = replace(String(output_path), r"\.html?$"i => "") * ".data.zz"
    sidecar_name = basename(sidecar_path)
    sidecar_name_js = replace(sidecar_name, "\\" => "\\\\", "\"" => "\\\"")
    raw_bytes = Vector{UInt8}(codeunits(runs_master_json))
    compressed_bytes = Zarr.zcompress(raw_bytes, Zarr.ZlibCompressor(7))
    write(sidecar_path, compressed_bytes)

    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>$(title)</title>
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" integrity="sha256-p4NxAoJBhIIN+hmNHrzRCf9tD/miZyoHS5obTRR9BMY=" crossorigin=""/>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <style>
        :root {
            --bg-dark: #0b0f19;
            --panel-bg: rgba(15, 23, 42, 0.90);
            --panel-border: rgba(255, 255, 255, 0.12);
            --accent-blue: #38bdf8;
            --accent-cyan: #06b6d4;
            --accent-green: #10b981;
            --accent-orange: #f97316;
            --accent-red: #ef4444;
            --accent-purple: #a855f7;
            --accent-amber: #f59e0b;
            --text-main: #f1f5f9;
            --text-muted: #94a3b8;
            --font-sans: 'Inter', system-ui, -apple-system, sans-serif;
            --font-mono: 'JetBrains Mono', monospace;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }

        body {
            font-family: var(--font-sans);
            background: var(--bg-dark);
            color: var(--text-main);
            height: 100vh;
            display: flex;
            flex-direction: column;
            overflow: hidden;
        }

        header {
            background: rgba(11, 15, 25, 0.96);
            backdrop-filter: blur(12px);
            border-bottom: 1px solid var(--panel-border);
            padding: 10px 20px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            z-index: 1000;
        }

        .header-title h1 {
            font-size: 1.1rem;
            font-weight: 700;
            color: #fff;
            letter-spacing: -0.02em;
        }

        .header-title h1 span {
            background: linear-gradient(135deg, var(--accent-blue), #818cf8);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }

        .header-subtitle {
            font-size: 0.76rem;
            color: var(--text-muted);
            margin-top: 2px;
        }

        .header-controls {
            display: flex;
            align-items: center;
            gap: 16px;
        }

        .run-selector-box {
            display: flex;
            align-items: center;
            gap: 8px;
            background: rgba(255, 255, 255, 0.06);
            padding: 4px 10px;
            border-radius: 8px;
            border: 1px solid var(--panel-border);
        }

        .run-selector-box label {
            font-size: 0.72rem;
            text-transform: uppercase;
            font-weight: 600;
            color: var(--accent-blue);
        }

        .run-selector-box select {
            background: #0f172a;
            color: #fff;
            border: 1px solid rgba(255,255,255,0.2);
            border-radius: 6px;
            padding: 4px 8px;
            font-size: 0.82rem;
            font-weight: 500;
            cursor: pointer;
            outline: none;
        }

        .metric-chips {
            display: flex;
            gap: 10px;
            align-items: center;
        }

        .chip {
            background: var(--panel-bg);
            border: 1px solid var(--panel-border);
            border-radius: 8px;
            padding: 4px 10px;
            display: flex;
            flex-direction: column;
            align-items: flex-end;
        }

        .chip-label {
            font-size: 0.62rem;
            text-transform: uppercase;
            letter-spacing: 0.05em;
            color: var(--text-muted);
        }

        .chip-value {
            font-size: 0.88rem;
            font-weight: 600;
            font-family: var(--font-mono);
            color: var(--accent-blue);
        }

        .chip-value.success { color: var(--accent-green); }
        .chip-value.warning { color: var(--accent-orange); }

        #workspace {
            position: relative;
            flex: 1;
            display: flex;
            width: 100%;
            height: calc(100vh - 60px);
        }

        #map {
            width: 100%;
            height: 100%;
            background: #080d1a;
        }

        .floating-panel {
            position: absolute;
            background: var(--panel-bg);
            backdrop-filter: blur(16px);
            border: 1px solid var(--panel-border);
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.45);
            z-index: 1000;
            padding: 14px;
            transition: all 0.2s ease;
        }

        #layer-manager-dock {
            top: 20px;
            right: 20px;
            width: 270px;
            max-height: calc(100vh - 200px);
            overflow-y: auto;
        }

        #layer-manager-dock::-webkit-scrollbar { width: 4px; }
        #layer-manager-dock::-webkit-scrollbar-thumb { background: rgba(255,255,255,0.2); border-radius: 2px; }

        .layer-group-title {
            font-size: 0.72rem;
            text-transform: uppercase;
            color: var(--accent-blue);
            margin: 10px 0 6px 0;
            letter-spacing: 0.06em;
            font-weight: 700;
            border-bottom: 1px solid rgba(255,255,255,0.08);
            padding-bottom: 3px;
        }

        .layer-group-title:first-child { margin-top: 0; }

        .layer-toggle-group {
            display: flex;
            flex-direction: column;
            gap: 6px;
        }

        .layer-item {
            display: flex;
            align-items: center;
            justify-content: space-between;
            font-size: 0.78rem;
            color: var(--text-main);
            padding: 3px 0;
        }

        .layer-item label {
            display: flex;
            align-items: center;
            gap: 6px;
            cursor: pointer;
            user-select: none;
        }

        .layer-item input[type="checkbox"] {
            accent-color: var(--accent-blue);
            cursor: pointer;
        }

        #controls-dock {
            top: 20px;
            left: 20px;
            width: 260px;
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .control-section h4 {
            font-size: 0.70rem;
            text-transform: uppercase;
            letter-spacing: 0.06em;
            color: var(--accent-blue);
            margin-bottom: 5px;
        }

        .control-select {
            width: 100%;
            background: rgba(0, 0, 0, 0.4);
            border: 1px solid var(--panel-border);
            color: #fff;
            padding: 5px 8px;
            border-radius: 6px;
            font-size: 0.78rem;
            outline: none;
        }

        .slider-control-row {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 8px;
            font-size: 0.74rem;
            color: var(--text-muted);
        }

        #playback-dock {
            bottom: 24px;
            left: 50%;
            transform: translateX(-50%);
            width: 90%;
            max-width: 720px;
            display: flex;
            flex-direction: column;
            gap: 10px;
        }

        .playback-controls {
            display: flex;
            align-items: center;
            justify-content: space-between;
            gap: 12px;
        }

        .btn-group { display: flex; gap: 6px; }

        .control-btn {
            background: rgba(255, 255, 255, 0.08);
            border: 1px solid var(--panel-border);
            color: #fff;
            padding: 6px 14px;
            border-radius: 6px;
            font-size: 0.82rem;
            font-weight: 500;
            cursor: pointer;
            display: inline-flex;
            align-items: center;
            gap: 6px;
            transition: background 0.15s ease;
        }

        .control-btn:hover { background: rgba(255, 255, 255, 0.18); }
        .control-btn.active {
            background: var(--accent-blue);
            color: #000;
            font-weight: 600;
            border-color: var(--accent-blue);
        }

        .time-display {
            font-family: var(--font-mono);
            font-size: 0.88rem;
            color: #fff;
            background: rgba(0, 0, 0, 0.35);
            padding: 4px 10px;
            border-radius: 6px;
            border: 1px solid var(--panel-border);
        }

        .slider-row { display: flex; align-items: center; gap: 12px; }

        input[type="range"] {
            flex: 1;
            accent-color: var(--accent-blue);
            height: 6px;
            border-radius: 3px;
            background: rgba(255, 255, 255, 0.2);
            cursor: pointer;
        }

        #hud-dock {
            bottom: 24px;
            right: 20px;
            width: 260px;
            font-size: 0.78rem;
        }

        .hud-grid {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 6px 10px;
            margin-top: 6px;
        }

        .hud-item { display: flex; flex-direction: column; }
        .hud-label { font-size: 0.62rem; color: var(--text-muted); text-transform: uppercase; }
        .hud-val { font-family: var(--font-mono); font-weight: 600; color: #fff; }

        #legend-dock {
            bottom: 125px;
            left: 20px;
            width: 250px;
            padding: 10px 12px;
            font-size: 0.74rem;
        }

        .legend-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-weight: 600;
            color: var(--accent-blue);
            margin-bottom: 6px;
            text-transform: uppercase;
            font-size: 0.68rem;
            letter-spacing: 0.05em;
        }

        .legend-bar {
            height: 12px;
            border-radius: 4px;
            width: 100%;
            margin-bottom: 4px;
            border: 1px solid rgba(255,255,255,0.2);
        }

        .legend-labels {
            display: flex;
            justify-content: space-between;
            font-family: var(--font-mono);
            font-size: 0.68rem;
            color: var(--text-muted);
        }

        .leaflet-tooltip {
            background: rgba(15, 23, 42, 0.94) !important;
            color: #fff !important;
            border: 1px solid rgba(255,255,255,0.25) !important;
            border-radius: 6px !important;
            font-family: var(--font-sans) !important;
            font-size: 0.76rem !important;
            box-shadow: 0 4px 16px rgba(0,0,0,0.5) !important;
            padding: 6px 10px !important;
            line-height: 1.35 !important;
        }

        #data-loader-overlay {
            position: fixed;
            inset: 0;
            z-index: 2000;
            display: flex;
            align-items: center;
            justify-content: center;
            background: rgba(11, 15, 25, 0.86);
            backdrop-filter: blur(6px);
        }

        #data-loader-overlay.hidden { display: none; }

        .loader-card {
            background: var(--panel-bg);
            border: 1px solid var(--panel-border);
            border-radius: 10px;
            padding: 20px 24px;
            max-width: 440px;
            width: 90%;
            box-shadow: 0 12px 32px rgba(0, 0, 0, 0.6);
            text-align: center;
        }

        .loader-title {
            font-size: 0.95rem;
            font-weight: 600;
            color: var(--accent-blue);
            margin-bottom: 8px;
        }

        .loader-status {
            font-size: 0.82rem;
            color: var(--text-muted);
            line-height: 1.45;
            margin-bottom: 12px;
        }

        .loader-btn {
            display: inline-block;
            padding: 8px 14px;
            border-radius: 6px;
            border: 1px solid var(--accent-blue);
            background: rgba(56, 189, 248, 0.14);
            color: var(--text-main);
            font-family: var(--font-sans);
            font-size: 0.82rem;
            font-weight: 500;
            cursor: pointer;
        }

        .loader-btn:hover { background: rgba(56, 189, 248, 0.24); }
    </style>
</head>
<body>
    <header>
        <div class="header-title">
            <h1>Larval <span>Dispersal</span> & Demographic Connectivity</h1>
            <div class="header-subtitle">Multi-Run Lagrangian Particle Tracking Dashboard (*Chionoecetes opilio*)</div>
        </div>
        <div class="header-controls">
            <div class="run-selector-box">
                <label for="run-select">Run:</label>
                <select id="run-select" onchange="switchRun(this.value)"></select>
            </div>
            <div class="metric-chips">
                <div class="chip">
                    <span class="chip-label">Cohort</span>
                    <span class="chip-value" id="chip-cohort">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Settled</span>
                    <span class="chip-value success" id="chip-settled">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Survival</span>
                    <span class="chip-value" id="chip-survival">-</span>
                </div>
                <div class="chip">
                    <span class="chip-label">Mean PLD</span>
                    <span class="chip-value warning" id="chip-pld">-</span>
                </div>
            </div>
        </div>
    </header>

    <div id="workspace">
        <div id="map"></div>

        <div class="floating-panel" id="controls-dock">
            <div class="control-section">
                <h4>Particle Palette</h4>
                <select id="color-mode-select" class="control-select">
                    <option value="stage">Larval Stage (Molting)</option>
                    <option value="depth">Drift Depth (m)</option>
                    <option value="temp">Ambient Temperature (°C)</option>
                    <option value="dd">Cumulative Degree-Days</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Larval Cohort Filter</h4>
                <select id="outcome-filter-select" class="control-select">
                    <option value="all">All Particles</option>
                    <option value="settled">Settled Successfully</option>
                    <option value="pelagic">Remaining Pelagic</option>
                    <option value="dead">Thermal Mortality (Dead)</option>
                </select>
            </div>

            <div class="control-section">
                <h4>Focus Particle</h4>
                <select id="particle-select" class="control-select">
                    <option value="-1">Overview (All Larvae)</option>
                </select>
            </div>
        </div>

        <div class="floating-panel" id="layer-manager-dock">
            <div class="layer-group-title">🛰️ Lagrangian Tracking</div>
            <div class="layer-toggle-group">
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-fulltracks" checked> 🛰️ Full Trajectories</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-progtrails" checked> ✨ Progressive Trails</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-heads" checked> 🔵 Active Larval Heads</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-release" checked> 🟢 Spawning Releases</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-settle" checked> 🔴 Settlement Sites</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-density"> 🔥 Nursery Settlement Density</label>
                </div>
            </div>

            <div class="layer-group-title">🏛️ Demographics & Zones</div>
            <div class="layer-toggle-group">
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-strata" checked> 🗺️ CFA Management Zones</label>
                </div>
                <div class="layer-item">
                    <label><input type="checkbox" id="layer-conn"> ↔️ CFA Connectivity Flows</label>
                </div>
            </div>
        </div>

        <div class="floating-panel" id="playback-dock">
            <div class="playback-controls">
                <div class="btn-group">
                    <button class="control-btn" id="play-btn">▶ Play</button>
                    <button class="control-btn" id="reset-btn">⏮ Reset</button>
                    <button class="control-btn speed-btn active" data-speed="1">1x</button>
                    <button class="control-btn speed-btn" data-speed="2">2x</button>
                    <button class="control-btn speed-btn" data-speed="4">4x</button>
                </div>
                <div class="time-display" id="time-readout">Day 0.00 / 0.00 d</div>
            </div>
            <div class="slider-row">
                <input type="range" id="time-slider" min="0" max="0" value="0" step="1">
            </div>
        </div>

        <div class="floating-panel" id="hud-dock">
            <div style="display:flex; justify-content:space-between; align-items:center;">
                <b style="color:var(--accent-blue);" id="telemetry-hud">Telemetry HUD</b>
                <span id="hud-id" style="font-family:var(--font-mono); color:#fff;">Overview</span>
            </div>
            <div class="hud-grid">
                <div class="hud-item"><span class="hud-label">Stage</span><span class="hud-val" id="hud-stage">-</span></div>
                <div class="hud-item"><span class="hud-label">Depth</span><span class="hud-val" id="hud-depth">-</span></div>
                <div class="hud-item"><span class="hud-label">Temp (T)</span><span class="hud-val" id="hud-temp">-</span></div>
                <div class="hud-item"><span class="hud-label">Degree-Days</span><span class="hud-val" id="hud-dd">-</span></div>
                <div class="hud-item"><span class="hud-label">Survival</span><span class="hud-val" id="hud-surv">-</span></div>
                <div class="hud-item"><span class="hud-label">Position</span><span class="hud-val" id="hud-pos">-</span></div>
            </div>
        </div>

        <div class="floating-panel" id="legend-dock">
            <div class="legend-header">
                <span id="legend-title">Nursery Settlement Density</span>
                <span id="legend-unit" style="font-family:var(--font-mono)">%</span>
            </div>
            <div class="legend-bar" id="legend-bar" style="background: linear-gradient(to right, #1e3a8a, #06b6d4, #f59e0b, #ef4444);"></div>
            <div class="legend-labels">
                <span id="legend-min">0.0</span>
                <span id="legend-mid">50.0</span>
                <span id="legend-max">100.0</span>
            </div>
        </div>
    </div>

    <div id="data-loader-overlay">
        <div class="loader-card">
            <div class="loader-title">Loading Trajectory Payload</div>
            <div class="loader-status" id="loader-status">Fetching $(sidecar_name)...</div>
            <div id="loader-picker" style="display:none;">
                <input type="file" id="sidecar-file-input" accept=".zz,.bin,.gz" style="display:none;" />
                <label for="sidecar-file-input" class="loader-btn">Select $(sidecar_name)</label>
            </div>
        </div>
    </div>

    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js" integrity="sha256-20nQCchB9co0qIjJZRGuk2/Z9VM+kNiyxNV1lvTlZBo=" crossorigin=""></script>
    <script>
        const SIDECAR_FILE = "$(sidecar_name_js)";
        const STRATA = $(strata_json);
        let RUNS = {};
        let activeRunKey = "$(default_run_key)";
        let currentRun = null;

        const STAGE_COLORS = {
            "zoea1": "#38bdf8", "zoea2": "#06b6d4", "zoea3": "#0d9488", "zoea4": "#eab308",
            "megalopa": "#f97316", "instar1_settled": "#10b981", "dead": "#ef4444"
        };

        const map = L.map('map', { center: [44.2, -62.5], zoom: 7 });
        const esriOcean = L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}', { maxZoom: 13 }).addTo(map);
        const cartoDark = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', { maxZoom: 19 });

        L.control.layers({ "ESRI Ocean": esriOcean, "CartoDB Dark": cartoDark }, null, { position: 'bottomleft' }).addTo(map);

        // Management & Lagrangian Layers
        const strataLayer = L.layerGroup().addTo(map);
        const fullTracksLayer = L.layerGroup().addTo(map);
        const progressTracksLayer = L.layerGroup().addTo(map);
        const activeHeadsLayer = L.layerGroup().addTo(map);
        const releaseMarkersLayer = L.layerGroup().addTo(map);
        const settlementMarkersLayer = L.layerGroup().addTo(map);
        const settlementDensityLayer = L.layerGroup();
        const connectivityFlowsLayer = L.layerGroup();

        STRATA.forEach(s => {
            let layer = s.is_polygon ? L.polygon(s.polygon, {color: s.color, weight: 2.2, fillOpacity: 0.08}) : L.rectangle([[s.lat_min, s.lon_min], [s.lat_max, s.lon_max]], {color: s.color, weight: 2, fillOpacity: 0.06});
            layer.bindTooltip("<b>" + s.name + "</b>").addTo(strataLayer);
        });

        const runSelect = document.getElementById('run-select');

        let currentStep = 0, isPlaying = false, playSpeed = 1, animTimer = null, showTrails = true, colorMode = "stage", outcomeFilter = "all", focusParticleId = -1;
        let headMarkers = [], fullTrackPolylines = [], progressTrackPolylines = [];

        function getColor(p, step) {
            if (colorMode === "stage") return STAGE_COLORS[p.stages[step]] || "#38bdf8";
            if (colorMode === "depth") return Math.abs(p.depths[step]) < 50 ? "#38bdf8" : "#1e3a8a";
            if (colorMode === "temp") return p.temps[step] < 6 ? "#38bdf8" : "#f87171";
            return p.dds[step] < 100 ? "#38bdf8" : "#f97316";
        }

        function updateLegend(title, unit, minV, maxV, gradientCss) {
            document.getElementById('legend-title').textContent = title;
            document.getElementById('legend-unit').textContent = unit;
            document.getElementById('legend-min').textContent = minV;
            document.getElementById('legend-max').textContent = maxV;
            document.getElementById('legend-mid').textContent = ((parseFloat(minV) + parseFloat(maxV)) / 2.0).toFixed(1);
            document.getElementById('legend-bar').style.background = gradientCss;
        }

        function renderActiveRun() {
            [fullTracksLayer, progressTracksLayer, activeHeadsLayer, releaseMarkersLayer,
             settlementMarkersLayer, settlementDensityLayer, connectivityFlowsLayer].forEach(l => l.clearLayers());

            headMarkers = [];
            fullTrackPolylines = [];
            progressTrackPolylines = [];

            const pList = currentRun.particles;
            const nTimes = currentRun.times.length;

            // 1. Summary Chips
            document.getElementById('chip-cohort').textContent = currentRun.stats.n_particles;
            document.getElementById('chip-settled').textContent = `\${currentRun.stats.n_settled} (\${currentRun.stats.settle_pct}%)`;
            document.getElementById('chip-survival').textContent = `\${currentRun.stats.survival_pct}%`;
            document.getElementById('chip-pld').textContent = `\${currentRun.stats.mean_pld_days} d`;

            // 2. Trajectories, Heads & Markers
            pList.forEach((p) => {
                const fullLatLngs = [];
                for (let t = 0; t < nTimes; t++) {
                    fullLatLngs.push([p.lats[t], p.lons[t]]);
                }

                // Full Polyline
                const fullPoly = L.polyline(fullLatLngs, {
                    color: getColor(p, nTimes - 1),
                    weight: 2.2,
                    opacity: 0.55
                }).on('click', () => selectParticle(p.id));
                fullTrackPolylines.push(fullPoly);
                fullTracksLayer.addLayer(fullPoly);

                // Progress Polyline
                const progPoly = L.polyline(fullLatLngs.slice(0, 2), {
                    color: getColor(p, 0),
                    weight: 3.5,
                    opacity: 0.95
                }).on('click', () => selectParticle(p.id));
                progressTrackPolylines.push(progPoly);
                progressTracksLayer.addLayer(progPoly);

                // Head Marker
                const marker = L.circleMarker([p.lats[0], p.lons[0]], {
                    radius: 5.5,
                    color: '#ffffff',
                    weight: 2,
                    fillColor: getColor(p, 0),
                    fillOpacity: 1.0
                }).on('click', () => selectParticle(p.id));
                marker.bindTooltip(`<b>Larva #\${p.id}</b> (\${p.stages[0]})`, { sticky: true });
                headMarkers.push(marker);
                activeHeadsLayer.addLayer(marker);

                // Spawning Release Marker
                const relMarker = L.circleMarker([p.lats[0], p.lons[0]], {
                    radius: 4,
                    color: '#10b981',
                    fillColor: '#059669',
                    fillOpacity: 0.85,
                    weight: 1.5
                }).bindPopup(`<b>Spawning Release #\${p.id}</b><br>Lon: \${p.lons[0].toFixed(3)}°E<br>Lat: \${p.lats[0].toFixed(3)}°N`);
                releaseMarkersLayer.addLayer(relMarker);

                // Settlement Marker (if settled)
                if (p.status === "settled_successful") {
                    const lastT = nTimes - 1;
                    const setMarker = L.circleMarker([p.lats[lastT], p.lons[lastT]], {
                        radius: 6,
                        color: '#f59e0b',
                        fillColor: '#d97706',
                        fillOpacity: 0.9,
                        weight: 2
                    }).bindPopup(`<b>Benthic Settlement Site #\${p.id}</b><br>Age: \${p.settlement_age_days} d<br>Lon: \${p.lons[lastT].toFixed(3)}°E<br>Lat: \${p.lats[lastT].toFixed(3)}°N`);
                    settlementMarkersLayer.addLayer(setMarker);
                }
            });

            // 3. Nursery Settlement Density Heatmap
            if (currentRun.settlement_density && currentRun.settlement_density.length > 0) {
                const maxDens = Math.max(...currentRun.settlement_density.map(d => d.density));
                currentRun.settlement_density.forEach(d => {
                    const norm = maxDens > 0 ? (d.density / maxDens) : 0.5;
                    const col = norm > 0.66 ? '#ef4444' : norm > 0.33 ? '#f59e0b' : '#38bdf8';
                    const rect = L.circle([d.lat, d.lon], {
                        radius: 7000,
                        color: col,
                        fillColor: col,
                        fillOpacity: 0.35 + norm * 0.45,
                        weight: 1
                    }).bindTooltip(`<b>Settlement Nursery Density</b><br>Probability: \${(d.density * 100).toFixed(2)}%`);
                    settlementDensityLayer.addLayer(rect);
                });
            }

            // 5. Demographic Connectivity Flows
            if (currentRun.connectivity_flows && currentRun.connectivity_flows.length > 0) {
                currentRun.connectivity_flows.forEach(flow => {
                    const weight = Math.max(2.0, flow.prob * 10.0);
                    const flowLine = L.polyline([[flow.src_lat, flow.src_lon], [flow.dst_lat, flow.dst_lon]], {
                        color: '#a855f7',
                        weight: weight,
                        opacity: 0.8,
                        dashArray: '8, 4'
                    }).bindTooltip(`<b>Connectivity Transition</b><br>\${flow.src} ➔ \${flow.dst}<br>Probability P_{ij}: \${(flow.prob * 100).toFixed(1)}%`);
                    connectivityFlowsLayer.addLayer(flowLine);
                });
            }

            // 6. Particle Selector Dropdown
            const pSelect = document.getElementById('particle-select');
            pSelect.innerHTML = '<option value="-1">Overview (All Larvae)</option>';
            pList.forEach(p => {
                const opt = document.createElement('option');
                opt.value = p.id;
                opt.textContent = `Larva #\${p.id} (\${p.status})`;
                pSelect.appendChild(opt);
            });

            // 7. Slider Bounds
            const slider = document.getElementById('time-slider');
            slider.max = nTimes - 1;
            slider.value = 0;

            // Auto-fit Bounds (iterative min/max to avoid V8 call-stack limits on large cohorts)
            let minLat = Infinity, maxLat = -Infinity, minLon = Infinity, maxLon = -Infinity;
            for (const p of pList) {
                for (let t = 0; t < p.lats.length; t++) {
                    const la = p.lats[t], lo = p.lons[t];
                    if (la !== null && Number.isFinite(la) && lo !== null && Number.isFinite(lo)) {
                        if (la < minLat) minLat = la;
                        if (la > maxLat) maxLat = la;
                        if (lo < minLon) minLon = lo;
                        if (lo > maxLon) maxLon = lo;
                    }
                }
            }
            if (Number.isFinite(minLat)) {
                map.fitBounds([
                    [minLat - 0.25, minLon - 0.25],
                    [maxLat + 0.25, maxLon + 0.25]
                ]);
            }

            updateLegend("Nursery Settlement Density", "%", "0.0", "100.0", "linear-gradient(to right, #1e3a8a, #06b6d4, #f59e0b, #ef4444)");
            updateFrame(0);
        }

        function switchRun(key) {
            if (RUNS[key]) {
                activeRunKey = key;
                currentRun = RUNS[key];
                isPlaying = false;
                document.getElementById('play-btn').textContent = "▶ Play";
                clearTimeout(animTimer);
                renderActiveRun();
            }
        }

        function updateFrame(step) {
            currentStep = step;
            const nTimes = currentRun.times.length;
            const pList = currentRun.particles;

            document.getElementById('time-slider').value = step;
            document.getElementById('time-readout').textContent = `Day \${currentRun.times[step].toFixed(2)} / \${currentRun.times[nTimes - 1].toFixed(2)} d`;

            pList.forEach((p, idx) => {
                const isFiltered = (outcomeFilter === "settled" && p.status !== "settled_successful") ||
                                   (outcomeFilter === "pelagic" && p.status !== "pelagic") ||
                                   (outcomeFilter === "dead" && p.alive);
                const isFocused = focusParticleId === -1 || focusParticleId === p.id;

                if (isFiltered || !isFocused) {
                    if (headMarkers[idx]) headMarkers[idx].setStyle({ opacity: 0, fillOpacity: 0 });
                    if (fullTrackPolylines[idx]) fullTrackPolylines[idx].setStyle({ opacity: 0 });
                    if (progressTrackPolylines[idx]) progressTrackPolylines[idx].setStyle({ opacity: 0 });
                    return;
                }

                const lat = p.lats[step];
                const lon = p.lons[step];
                const col = getColor(p, step);

                if (headMarkers[idx]) {
                    headMarkers[idx].setLatLng([lat, lon]);
                    headMarkers[idx].setStyle({
                        fillColor: col,
                        opacity: 1,
                        fillOpacity: 1,
                        radius: (focusParticleId === p.id) ? 9 : 5.5
                    });
                    headMarkers[idx].setTooltipContent(`<b>Larva #\${p.id}</b> (\${p.stages[step]} • \${p.depths[step].toFixed(1)}m)`);
                }

                if (fullTrackPolylines[idx]) {
                    fullTrackPolylines[idx].setStyle({
                        color: col,
                        opacity: (focusParticleId === p.id) ? 0.9 : 0.45,
                        weight: (focusParticleId === p.id) ? 3.5 : 2.0
                    });
                }

                if (progressTrackPolylines[idx]) {
                    if (showTrails && step > 0) {
                        const progPts = [];
                        for (let t = 0; t <= step; t++) {
                            progPts.push([p.lats[t], p.lons[t]]);
                        }
                        progressTrackPolylines[idx].setLatLngs(progPts);
                        progressTrackPolylines[idx].setStyle({
                            color: col,
                            opacity: (focusParticleId === p.id) ? 1.0 : 0.95,
                            weight: (focusParticleId === p.id) ? 4.5 : 3.2
                        });
                    } else if (showTrails && step === 0) {
                        const initPts = [[p.lats[0], p.lons[0]], [p.lats[Math.min(1, nTimes - 1)], p.lons[Math.min(1, nTimes - 1)]]];
                        progressTrackPolylines[idx].setLatLngs(initPts);
                        progressTrackPolylines[idx].setStyle({ color: col, opacity: 0.9, weight: 3.0 });
                    } else {
                        progressTrackPolylines[idx].setStyle({ opacity: 0 });
                    }
                }
            });

            if (focusParticleId !== -1) {
                const p = pList.find(item => item.id === focusParticleId);
                if (p) updateHUD(p, step);
            }
        }

        function updateHUD(p, step) {
            document.getElementById('hud-id').textContent = `#\${p.id}`;
            document.getElementById('hud-stage').textContent = p.stages[step];
            document.getElementById('hud-depth').textContent = p.depths[step] !== null ? `\${p.depths[step].toFixed(1)} m` : '-';
            document.getElementById('hud-temp').textContent = p.temps[step] !== null ? `\${p.temps[step].toFixed(2)} °C` : '-';
            document.getElementById('hud-dd').textContent = p.dds[step] !== null ? `\${p.dds[step].toFixed(1)} DD` : '-';
            document.getElementById('hud-surv').textContent = p.survs[step] !== null ? `\${(p.survs[step] * 100).toFixed(1)} %` : '-';
            document.getElementById('hud-pos').textContent =
                (p.lats[step] !== null && p.lons[step] !== null)
                    ? `\${p.lats[step].toFixed(2)}°N, \${p.lons[step].toFixed(2)}°E`
                    : '-';
        }

        function selectParticle(id) {
            focusParticleId = id;
            document.getElementById('particle-select').value = id;
            const p = currentRun.particles.find(item => item.id === id);
            if (p) {
                map.setView([p.lats[currentStep], p.lons[currentStep]], 9);
                updateHUD(p, currentStep);
            }
            updateFrame(currentStep);
        }

        function playLoop() {
            if (!isPlaying) return;
            const nTimes = currentRun.times.length;
            let nextStep = currentStep + playSpeed;
            if (nextStep >= nTimes) nextStep = 0;
            updateFrame(nextStep);
            animTimer = setTimeout(playLoop, 60);
        }

        // Cursor Coordinate Readout
        map.on('mousemove', function(e) {
            if (focusParticleId !== -1) return;
            document.getElementById('hud-id').textContent =
                `\${e.latlng.lat.toFixed(2)}°N, \${e.latlng.lng.toFixed(2)}°E`;
        });

        // Layer Manager Event Listeners
        const layerControls = [
            { id: 'layer-fulltracks', layer: fullTracksLayer },
            { id: 'layer-progtrails', layer: progressTracksLayer },
            { id: 'layer-heads', layer: activeHeadsLayer },
            { id: 'layer-release', layer: releaseMarkersLayer },
            { id: 'layer-settle', layer: settlementMarkersLayer },
            { id: 'layer-density', layer: settlementDensityLayer, onSelect: () => updateLegend("Nursery Settlement Density", "%", "0.0", "100.0", "linear-gradient(to right, #1e3a8a, #06b6d4, #f59e0b, #ef4444)") },
            { id: 'layer-strata', layer: strataLayer },
            { id: 'layer-conn', layer: connectivityFlowsLayer }
        ];

        layerControls.forEach(ctrl => {
            const el = document.getElementById(ctrl.id);
            if (el) {
                el.addEventListener('change', function() {
                    if (this.checked) {
                        map.addLayer(ctrl.layer);
                        if (ctrl.onSelect) ctrl.onSelect();
                    } else {
                        map.removeLayer(ctrl.layer);
                    }
                });
            }
        });

        // Controller Event Listeners
        document.getElementById('play-btn').addEventListener('click', function() {
            isPlaying = !isPlaying;
            this.textContent = isPlaying ? "⏸ Pause" : "▶ Play";
            if (isPlaying) playLoop();
            else clearTimeout(animTimer);
        });

        document.getElementById('reset-btn').addEventListener('click', function() {
            isPlaying = false;
            document.getElementById('play-btn').textContent = "▶ Play";
            clearTimeout(animTimer);
            updateFrame(0);
        });

        document.getElementById('time-slider').addEventListener('input', function(e) {
            updateFrame(parseInt(e.target.value));
        });

        document.querySelectorAll('.speed-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                document.querySelectorAll('.speed-btn').forEach(b => b.classList.remove('active'));
                this.classList.add('active');
                playSpeed = parseInt(this.getAttribute('data-speed'));
            });
        });

        document.getElementById('color-mode-select').addEventListener('change', function(e) {
            colorMode = e.target.value;
            updateFrame(currentStep);
        });

        document.getElementById('outcome-filter-select').addEventListener('change', function(e) {
            outcomeFilter = e.target.value;
            updateFrame(currentStep);
        });

        document.getElementById('particle-select').addEventListener('change', function(e) {
            selectParticle(parseInt(e.target.value));
        });

        async function decodeZlibJson(stream) {
            if (typeof DecompressionStream === "undefined") {
                throw new Error("Browser does not support DecompressionStream.");
            }
            const decompressed = stream.pipeThrough(new DecompressionStream("deflate"));
            const text = await new Response(decompressed).text();
            return JSON.parse(text);
        }

        function bootWithRuns(runsObj) {
            RUNS = runsObj;
            const keys = Object.keys(RUNS);
            if (keys.length === 0) {
                throw new Error("Sidecar payload contains no simulation runs.");
            }
            if (!RUNS[activeRunKey]) {
                activeRunKey = keys[0];
            }
            currentRun = RUNS[activeRunKey];
            runSelect.innerHTML = "";
            keys.forEach(k => {
                const opt = document.createElement('option');
                opt.value = k;
                opt.textContent = k;
                if (k === activeRunKey) opt.selected = true;
                runSelect.appendChild(opt);
            });
            document.getElementById('data-loader-overlay').classList.add('hidden');
            renderActiveRun();
        }

        async function loadSidecarAndBoot() {
            const statusEl = document.getElementById('loader-status');
            const pickerEl = document.getElementById('loader-picker');
            const fileInput = document.getElementById('sidecar-file-input');

            fileInput.addEventListener('change', async function(e) {
                const file = e.target.files && e.target.files[0];
                if (!file) return;
                statusEl.textContent = `Decompressing \${file.name}...`;
                try {
                    const runsObj = await decodeZlibJson(file.stream());
                    bootWithRuns(runsObj);
                } catch (err) {
                    statusEl.textContent = `Failed to decode \${file.name}: \${err.message}`;
                }
            });

            try {
                const resp = await fetch(`./\${SIDECAR_FILE}`);
                if (!resp.ok) {
                    throw new Error(`HTTP \${resp.status}`);
                }
                statusEl.textContent = `Decompressing \${SIDECAR_FILE}...`;
                const runsObj = await decodeZlibJson(resp.body);
                bootWithRuns(runsObj);
            } catch (err) {
                statusEl.textContent =
                    `Automatic fetch of \${SIDECAR_FILE} was blocked or unavailable (\${err.message}). ` +
                    `Serve this folder over HTTP/HTTPS (e.g. GitHub Pages or a local server), ` +
                    `or select \${SIDECAR_FILE} directly below.`;
                pickerEl.style.display = "block";
            }
        }

        loadSidecarAndBoot();
    </script>
</body>
</html>
"""

    write(output_path, html_content)
    raw_mb = round(length(raw_bytes) / 1.0e6; digits = 2)
    cmp_mb = round(length(compressed_bytes) / 1.0e6; digits = 2)
    println(
        "Exported multi-layer interactive dashboard to $(output_path) " *
        "(sidecar: $(sidecar_path), $(cmp_mb) MB compressed from $(raw_mb) MB)"
    )
    return output_path
end

"""
    export_interactive_tracks_html(
        trajectories::NamedTuple,
        output_path::AbstractString;
        kwargs...
    ) -> String

Multiple-dispatch convenience overload accepting `trajectories` as the primary argument,
forwarding to the file-path first method signature.

# Mathematical & Physical Context
Visualizes 4D Lagrangian trajectory vectors \$\\mathbf{x}_p(t) = (\\\\\\\\lambda_p(t), \\phi_p(t), z_p(t))\$
co-registered with Eulerian hydrodynamic fields (advection currents \$\\\\\\\\\\\\\\boldsymbol{u}_h\$,
temperature \$T\$, salinity \$S\$, stratification \$N^2\$, and turbulent diffusivity \$\\\\\\\\\\\\\\kappa_v\$)
across discrete depth levels and simulation epochs.

# Inputs
- `trajectories::NamedTuple`: Particle tracking output containing `lons`, `lats`, `depths`, `times`.
- `output_path::AbstractString`: Output path for the generated standalone HTML file.
- `kwargs...`: Additional keyword arguments forwarded to the file-path primary signature.

# Outputs
- `String`: Canonical file path to the generated HTML dashboard.
"""
function export_interactive_tracks_html(
    trajectories::NamedTuple,
    output_path::AbstractString;
    kwargs...
)::String
    return export_interactive_tracks_html(
        output_path;
        trajectories = trajectories,
        kwargs...
    )
end

"""
    plot_interactive_trajectories_map(
        trajectories::Union{Nothing, NamedTuple} = nothing;
        runs_data = nothing,
        hydrodynamics = nothing,
        gridded_dispersal = nothing,
        connectivity = nothing,
        recruitment_metrics = nothing,
        output_path::AbstractString = "outputs/interactive_larval_tracks.html",
        bathymetry = nothing,
        strata_definitions = nothing,
        title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
    ) -> String

Convenience wrapper to generate and export an interactive Leaflet HTML map of particle tracks,
supporting multi-run data dictionaries, hydrodynamic Eulerian flow fields, and full spatial layer suites.
"""
function plot_interactive_trajectories_map(
    trajectories::Union{Nothing, NamedTuple} = nothing;
    runs_data = nothing,
    hydrodynamics = nothing,
    gridded_dispersal = nothing,
    connectivity = nothing,
    recruitment_metrics = nothing,
    output_path::AbstractString = "outputs/interactive_larval_tracks.html",
    bathymetry = nothing,
    strata_definitions = nothing,
    title::AbstractString = "Interactive Larval Dispersal & Demographic Connectivity Map"
)::String
    return export_interactive_tracks_html(
        output_path;
        trajectories = trajectories,
        runs_data = runs_data,
        hydrodynamics = hydrodynamics,
        gridded_dispersal = gridded_dispersal,
        connectivity = connectivity,
        recruitment_metrics = recruitment_metrics,
        bathymetry = bathymetry,
        strata_definitions = strata_definitions,
        title = title
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Hydrodynamic Field Video & GIF Animations (Oceananigans / CairoMakie)
# ─────────────────────────────────────────────────────────────────────────────

"""
    resolve_hydrodynamic_frame_data(
        hydro_frame::NamedTuple,
        variable::Symbol
    ) -> Tuple{Matrix{Float64}, Symbol, String}

Extract a 2D scalar spatial matrix, default colormap, and descriptive physical unit label
for a given hydrodynamic diagnostic variable at the selected discrete depth level.

# Inputs
- `hydro_frame::NamedTuple`: Hydrodynamic snapshot dataset extracted via
  `extract_hydrodynamic_dataset`.
- `variable::Symbol`: Target variable symbol (`:speed`, `:vorticity`, `:temperature`,
  `:salinity`, `:density`, `:stratification`, `:elevation`, `:w`).

# Outputs
- `Tuple{Matrix{Float64}, Symbol, String}`: 2D field matrix, CairoMakie colormap, label.
"""
function resolve_hydrodynamic_frame_data(
    hydro_frame::NamedTuple,
    variable::Symbol
)::Tuple{Matrix{Float64}, Symbol, String}
    k = hydro_frame.depth_index
    z_cur = (hasproperty(hydro_frame, :depths) && k <= length(hydro_frame.depths)) ?
        hydro_frame.depths[k] : -10.0
    bathy = hasproperty(hydro_frame, :bathymetry) ? hydro_frame.bathymetry : nothing

    mask_dry(mat::AbstractMatrix) = begin
        out = copy(mat)
        if !isnothing(bathy) && size(bathy) == size(out)
            # Immersed cells: bottom is shallower than current depth or dry land (bathy >= 0)
            is_dry = (bathy .>= 0.0) .| (z_cur .< bathy)
            out[is_dry] .= NaN
        elseif any(out .!= 0.0)
            out[out .== 0.0] .= NaN
        end
        return out
    end

    if variable in (:salinity, :S, :sal)
        return (mask_dry(Float64.(hydro_frame.salinity[:, :, k])), :haline, "Practical Salinity (PSU)")
    elseif variable in (:temperature, :T, :temp, :theta)
        return (mask_dry(Float64.(hydro_frame.temperature[:, :, k])), :thermal, "Temperature (°C)")
    elseif variable in (:advection, :speed, :current, :velocity, :u_h)
        return (mask_dry(Float64.(hydro_frame.speed[:, :, k] .* 100.0)), :viridis, "Current Speed |u_h| (cm s⁻¹)")
    elseif variable in (:diffusion, :diffusivity, :kappa, :eddy_diffusion)
        return (mask_dry(Float64.(hydro_frame.diffusion[:, :, k] .* 10000.0)), :turbid, "Eddy Diffusivity κ (10⁻⁴ m² s⁻¹)")
    elseif variable in (:viscosity, :nu, :eddy_viscosity)
        return (mask_dry(Float64.(hydro_frame.viscosity[:, :, k] .* 10000.0)), :deep, "Eddy Viscosity ν (10⁻⁴ m² s⁻¹)")
    elseif variable in (:density, :rho)
        return (mask_dry(Float64.(hydro_frame.density[:, :, k])), :dense, "Potential Density (kg m⁻³)")
    elseif variable in (:stratification, :N2)
        return (mask_dry(Float64.(hydro_frame.stratification[:, :, k] .* 10000.0)), :ice, "Stratification N² (10⁻⁴ s⁻²)")
    elseif variable in (:richardson, :Ri)
        return (mask_dry(Float64.(hydro_frame.richardson_number[:, :, k])), :spectral, "Richardson Number Ri")
    elseif variable in (:vorticity, :zeta)
        vort_mat = ndims(hydro_frame.vorticity) == 2 ?
            Float64.(hydro_frame.vorticity[:, :]) :
            Float64.(hydro_frame.vorticity[:, :, k])
        return (mask_dry(vort_mat .* 100000.0), :balance, "Relative Vorticity ζ (10⁻⁵ s⁻¹)")
    elseif variable in (:elevation, :eta, :ssh)
        return (Float64.(hydro_frame.elevation .* 100.0), :delta, "Sea Surface Elevation η (cm)")
    elseif variable in (:w, :vertical_velocity)
        return (mask_dry(Float64.(hydro_frame.w[:, :, k] .* 1000.0)), :curl, "Vertical Velocity w (mm s⁻¹)")
    else
        return (mask_dry(Float64.(hydro_frame.speed[:, :, k] .* 100.0)), :viridis, "Current Speed |u_h| (cm s⁻¹)")
    end
end

