using Dates
using Random
using Statistics

"""
    DEFAULT_RUN_ID

The identifier a run is archived under when the configuration names none.

The store holds one run per configuration directory, so a fixed name is what says "this
run" rather than naming a climate pathway. A configuration that wants a second run beside
it names it in `[storage] run_id`.
"""
const DEFAULT_RUN_ID = "run"

"""
    run_identifier(opts::HydrodynamicOptions) -> String

The identifier this run is archived under in the trajectory store: the configured
`[storage] run_id` when there is one, and [`DEFAULT_RUN_ID`](@ref) otherwise.
"""
function run_identifier(opts::HydrodynamicOptions)::String
    tag = strip(opts.run_id)
    return isempty(tag) ? DEFAULT_RUN_ID : String(tag)
end

"""
    resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String) -> Tuple{String, String}

Resolve the target realized-simulation store path and name. If `opts.hydro_model_file` is
specified it is used directly (or resolved against the package root when it is a path,
or against `opts.output_dir` when it is a bare name); if empty it is
`joinpath(opts.output_dir, default_filename)`. The path is normally the
Tethys store the ocean half hands over as `--hydro-file`, though the name is left generic so
it also covers a JLD2 file a caller names.

# Inputs
- `opts::HydrodynamicOptions`: Configuration parameters.
- `default_filename::String`: Default fallback name when none is specified.

# Outputs
- `Tuple{String, String}`: `(full_path, file_basename)`
"""
function resolve_hydro_model_path(opts::HydrodynamicOptions, default_filename::String)
    if isempty(opts.hydro_model_file)
        full_path = joinpath(opts.output_dir, default_filename)
        return (full_path, default_filename)
    else
        raw = opts.hydro_model_file
        # A configured path is resolved against THIS PACKAGE'S ROOT, not against the
        # working directory the run happened to be started from. `../lakehouse/...` must
        # name the same store from every directory, and resolving it against the working
        # directory made the same configuration read a different simulation depending on
        # where it was run -- or none at all. A bare name is still this run's own output
        # directory, which is where a store this package produced lives.
        full_path = if isabspath(raw)
            raw
        elseif dirname(raw) != ""
            Tethys.project_lakehouse_path(raw, _PROJECT_ROOT)
        else
            joinpath(opts.output_dir, raw)
        end
        return (full_path, basename(raw))
    end
end

@inline function _bracket_axis(vals::Vector{Float64}, v::Float64, asc::Bool)
    n = length(vals)
    n <= 1 && return (1, 1, 0.0)
    idx = asc ? searchsortedfirst(vals, v) : searchsortedfirst(vals, v; rev = true)
    i = clamp(idx, 2, n)
    @inbounds v0 = vals[i - 1]
    @inbounds v1 = vals[i]
    denom = v1 - v0
    s = denom == 0.0 ? 0.0 : clamp((v - v0) / denom, 0.0, 1.0)
    return (i - 1, i, s)
end

@inline function _trilinear_f32(
    a::Array{Float32, 3},
    i0::Int, i1::Int, sx::Float64,
    j0::Int, j1::Int, sy::Float64,
    k0::Int, k1::Int, sz::Float64
)::Float64
    omx = 1.0 - sx
    omy = 1.0 - sy
    omz = 1.0 - sz
    w000 = omx * omy * omz
    w100 = sx  * omy * omz
    w010 = omx * sy  * omz
    w110 = sx  * sy  * omz
    w001 = omx * omy * sz
    w101 = sx  * omy * sz
    w011 = omx * sy  * sz
    w111 = sx  * sy  * sz
    @inbounds begin
        v000 = Float64(a[i0, j0, k0])
        v100 = Float64(a[i1, j0, k0])
        v010 = Float64(a[i0, j1, k0])
        v110 = Float64(a[i1, j1, k0])
        v001 = Float64(a[i0, j0, k1])
        v101 = Float64(a[i1, j0, k1])
        v011 = Float64(a[i0, j1, k1])
        v111 = Float64(a[i1, j1, k1])
    end
    res = w000 * v000 + w100 * v100 + w010 * v010 + w110 * v110 +
          w001 * v001 + w101 * v101 + w011 * v011 + w111 * v111
    isfinite(res) && return res

    # Near the seabed or coastline, deeper/land corners may be masked (NaN); renormalize
    # over the finite wet corners so near-bottom larvae sample the bottom-most wet layer.
    acc = 0.0
    wsum = 0.0
    if isfinite(v000); acc += w000 * v000; wsum += w000; end
    if isfinite(v100); acc += w100 * v100; wsum += w100; end
    if isfinite(v010); acc += w010 * v010; wsum += w010; end
    if isfinite(v110); acc += w110 * v110; wsum += w110; end
    if isfinite(v001); acc += w001 * v001; wsum += w001; end
    if isfinite(v101); acc += w101 * v101; wsum += w101; end
    if isfinite(v011); acc += w011 * v011; wsum += w011; end
    if isfinite(v111); acc += w111 * v111; wsum += w111; end
    return wsum > 0.0 ? acc / wsum : NaN
end

struct _FlowSlice
    u::Array{Float32, 3}
    v::Array{Float32, 3}
    w::Array{Float32, 3}
    T::Array{Float32, 3}
end

mutable struct _FlowWindowInterpolator
    store::HydrodynamicStore
    keep::Vector{Int}
    lon::Vector{Float64}
    lat::Vector{Float64}
    z::Vector{Float64}
    tsel::Vector{Float64}
    lon_asc::Bool
    lat_asc::Bool
    z_asc::Bool
    t_asc::Bool
    cached_m0::Int
    cached_m1::Int
    s0::_FlowSlice
    s1::_FlowSlice
end

function _read_flow_slice(store::HydrodynamicStore, idx::Int)::_FlowSlice
    return _FlowSlice(
        Float32.(variable(store, "u", idx)),
        Float32.(variable(store, "v", idx)),
        Float32.(variable(store, "w", idx)),
        Float32.(variable(store, "T", idx))
    )
end

function _ensure_flow_window!(fw::_FlowWindowInterpolator, m0::Int, m1::Int)
    if fw.cached_m0 == m0 && fw.cached_m1 == m1
        return nothing
    end
    if fw.cached_m1 == m0
        fw.s0 = fw.s1
        fw.cached_m0 = m0
        fw.s1 = (m1 == m0) ? fw.s0 : _read_flow_slice(fw.store, fw.keep[m1])
        fw.cached_m1 = m1
    else
        fw.s0 = _read_flow_slice(fw.store, fw.keep[m0])
        fw.cached_m0 = m0
        fw.s1 = (m1 == m0) ? fw.s0 : _read_flow_slice(fw.store, fw.keep[m1])
        fw.cached_m1 = m1
    end
    return nothing
end

function (fw::_FlowWindowInterpolator)(lon_q::Real, lat_q::Real, z_q::Real, t_q::Real)
    i0, i1, sx = _bracket_axis(fw.lon, Float64(lon_q), fw.lon_asc)
    j0, j1, sy = _bracket_axis(fw.lat, Float64(lat_q), fw.lat_asc)
    k0, k1, sz = _bracket_axis(fw.z, Float64(z_q), fw.z_asc)
    m0, m1, θ = _bracket_axis(fw.tsel, Float64(t_q), fw.t_asc)
    _ensure_flow_window!(fw, m0, m1)
    w0 = 1.0 - θ
    s0 = fw.s0
    s1 = fw.s1
    u = w0 * _trilinear_f32(s0.u, i0, i1, sx, j0, j1, sy, k0, k1, sz) +
        θ  * _trilinear_f32(s1.u, i0, i1, sx, j0, j1, sy, k0, k1, sz)
    v = w0 * _trilinear_f32(s0.v, i0, i1, sx, j0, j1, sy, k0, k1, sz) +
        θ  * _trilinear_f32(s1.v, i0, i1, sx, j0, j1, sy, k0, k1, sz)
    w = w0 * _trilinear_f32(s0.w, i0, i1, sx, j0, j1, sy, k0, k1, sz) +
        θ  * _trilinear_f32(s1.w, i0, i1, sx, j0, j1, sy, k0, k1, sz)
    T = w0 * _trilinear_f32(s0.T, i0, i1, sx, j0, j1, sy, k0, k1, sz) +
        θ  * _trilinear_f32(s1.T, i0, i1, sx, j0, j1, sy, k0, k1, sz)
    return (u = u, v = v, w = w, T = T)
end

function _sample_temperature(fw::_FlowWindowInterpolator, lon_q::Real, lat_q::Real, z_q::Real, t_q::Real)::Float64
    i0, i1, sx = _bracket_axis(fw.lon, Float64(lon_q), fw.lon_asc)
    j0, j1, sy = _bracket_axis(fw.lat, Float64(lat_q), fw.lat_asc)
    k0, k1, sz = _bracket_axis(fw.z, Float64(z_q), fw.z_asc)
    m0, m1, θ = _bracket_axis(fw.tsel, Float64(t_q), fw.t_asc)
    _ensure_flow_window!(fw, m0, m1)
    return (1.0 - θ) * _trilinear_f32(fw.s0.T, i0, i1, sx, j0, j1, sy, k0, k1, sz) +
           θ         * _trilinear_f32(fw.s1.T, i0, i1, sx, j0, j1, sy, k0, k1, sz)
end

"""
    _build_domain_interior_fn(fw::_FlowWindowInterpolator; margin_km = 25.0, rot_deg = 35.0)

Return a predicate `(lon, lat) -> Bool` that is `true` when `(lon, lat)` lies at least
`margin_km` inside the outer lateral edges of the realized hydrodynamic wet footprint
(measured in the shelf-aligned frame rotated `rot_deg` CCW from East). Used to stop
tracking larvae that approach open-ocean domain boundaries (`:off_domain`).
"""
function _build_domain_interior_fn(
    fw::_FlowWindowInterpolator;
    margin_km::Real = 25.0,
    rot_deg::Real = 35.0
)
    s0 = fw.s0
    nx, ny, _ = size(s0.u)
    k_surf = argmax(fw.z)
    lons = fw.lon
    lats = fw.lat

    lon_sum = 0.0
    lat_sum = 0.0
    wet_cnt = 0
    for j in 1:ny, i in 1:nx
        if isfinite(s0.u[i, j, k_surf]) && isfinite(s0.v[i, j, k_surf])
            lon_sum += lons[i]
            lat_sum += lats[j]
            wet_cnt += 1
        end
    end
    wet_cnt == 0 && return (_lo, _la) -> true
    lon0 = lon_sum / wet_cnt
    lat0 = lat_sum / wet_cnt

    km_per_lat = 111.132
    km_per_lon = 111.320 * cosd(lat0)
    c_r = cosd(Float64(rot_deg))
    s_r = sind(Float64(rot_deg))

    x_min, x_max = Inf, -Inf
    y_min, y_max = Inf, -Inf
    for j in 1:ny, i in 1:nx
        if isfinite(s0.u[i, j, k_surf]) && isfinite(s0.v[i, j, k_surf])
            dx = (lons[i] - lon0) * km_per_lon
            dy = (lats[j] - lat0) * km_per_lat
            xk = dx * c_r + dy * s_r
            yk = -dx * s_r + dy * c_r
            xk < x_min && (x_min = xk); xk > x_max && (x_max = xk)
            yk < y_min && (y_min = yk); yk > y_max && (y_max = yk)
        end
    end

    # Only enforce margins if the domain is meaningfully wider than 2 * margin_km.
    m = Float64(margin_km)
    if (x_max - x_min) <= 3.0 * m || (y_max - y_min) <= 3.0 * m
        return (_lo, _la) -> true
    end
    x_lo, x_hi = x_min + m, x_max - m
    y_lo, y_hi = y_min + m, y_max - m

    return (lo::Real, la::Real) -> begin
        dx = (Float64(lo) - lon0) * km_per_lon
        dy = (Float64(la) - lat0) * km_per_lat
        xk = dx * c_r + dy * s_r
        yk = -dx * s_r + dy * c_r
        return (x_lo <= xk <= x_hi) && (y_lo <= yk <= y_hi)
    end
end

function _build_static_wet_mask(
    fw::_FlowWindowInterpolator;
    bathy_elevation::Union{Nothing, AbstractMatrix} = nothing,
    coastline = nothing
)
    # Build a 2D horizontal wet-ocean mask from the shallowest (surface) model layer so
    # continental-shelf cells (where deeper z levels sit below the seabed) remain wet.
    s0 = fw.s0
    nx, ny, _ = size(s0.u)
    k_surf = argmax(fw.z)
    wet2d = BitMatrix(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        wet2d[i, j] = isfinite(s0.u[i, j, k_surf]) &&
                      isfinite(s0.v[i, j, k_surf]) &&
                      isfinite(s0.T[i, j, k_surf])
    end

    # Bake closed coastline polygons onto ALL wet grid cells once (regardless of smoothed
    # seabed depth, which can be < -50 m near steep coastal channels such as Cape North).
    if !isnothing(coastline) && !isempty(coastline)
        closed_polys = [
            (
                min_x = minimum(p.lons),
                max_x = maximum(p.lons),
                min_y = minimum(p.lats),
                max_y = maximum(p.lats),
                lons = p.lons,
                lats = p.lats
            ) for p in coastline if length(p.lons) >= 3 &&
                                    p.lons[begin] == p.lons[end] &&
                                    p.lats[begin] == p.lats[end]
        ]
        if !isempty(closed_polys)
            @inbounds for j in 1:ny
                la = fw.lat[j]
                for i in 1:nx
                    wet2d[i, j] || continue
                    lo = fw.lon[i]
                    for cp in closed_polys
                        if cp.min_x <= lo <= cp.max_x && cp.min_y <= la <= cp.max_y &&
                           Tethys.point_in_polygon(lo, la, cp.lons, cp.lats)
                            wet2d[i, j] = false
                            break
                        end
                    end
                end
            end
        end
    end

    lon_min, lon_max = extrema(fw.lon)
    lat_min, lat_max = extrema(fw.lat)
    lon = fw.lon
    lat = fw.lat
    lon_asc = fw.lon_asc
    lat_asc = fw.lat_asc

    return function wet_mask_fn(lon_q::Real, lat_q::Real, z_q::Real)::Bool
        lo = Float64(lon_q)
        la = Float64(lat_q)
        (lon_min <= lo <= lon_max && lat_min <= la <= lat_max) || return false
        i0, i1, _ = _bracket_axis(lon, lo, lon_asc)
        j0, j1, _ = _bracket_axis(lat, la, lat_asc)
        @inbounds return wet2d[i0, j0] && wet2d[i1, j0] && wet2d[i0, j1] && wet2d[i1, j1]
    end
end

"""
    build_flow_from_store(store; max_snapshots, start_time, span) -> _FlowWindowInterpolator

The advective field for tracking, read from an already-opened `HydrodynamicStore` via a
two-snapshot sliding window in memory. Only the two time slices bracketing `t_q` are held
as `Float32` 3D arrays at once (~340 MB for a 719 x 640 x 23 shelf grid rather than tens of
gigabytes for the full record), sliding forward without re-reading the shared endpoint as
`track_larval_cohort` advances monotonically in time. The retained span is written back
through `span`, and interpolation is trilinear in space and linear in time using O(log N)
binary search bracketing.
"""
function build_flow_from_store(store::HydrodynamicStore;
                               max_snapshots::Int = typemax(Int),
                               start_time::Real = -Inf,
                               span::Vector{Float64} = [0.0, 0.0])
    t = store.times
    valid_idx = findall(x -> x >= Float64(start_time), t)
    isempty(valid_idx) && error(
        "No hydrodynamic snapshots at or after start_time = $(start_time) s " *
        "(store spans $(t[1]) .. $(t[end]) s)."
    )
    nt = length(valid_idx)
    keep = (max_snapshots <= 0 || nt <= max_snapshots) ? valid_idx :
           valid_idx[sort(unique(round.(Int, range(1, stop = nt, length = max(2, max_snapshots)))))]
    tsel = Float64.(t[keep])

    lon = Float64.(collect(store.axes.lon))
    lat = Float64.(collect(store.axes.lat))
    z = Float64.(collect(store.axes.depth))

    s0 = _read_flow_slice(store, keep[1])
    m1_init = length(keep) > 1 ? 2 : 1
    s1 = m1_init == 1 ? s0 : _read_flow_slice(store, keep[m1_init])

    span[1] = tsel[1]
    span[2] = tsel[end]

    return _FlowWindowInterpolator(
        store,
        keep,
        lon,
        lat,
        z,
        tsel,
        length(lon) <= 1 || lon[end] >= lon[1],
        length(lat) <= 1 || lat[end] >= lat[1],
        length(z) <= 1 || z[end] >= z[1],
        length(tsel) <= 1 || tsel[end] >= tsel[1],
        1,
        m1_init,
        s0,
        s1
    )
end


"""
    run_segment_tracking(; opts::HydrodynamicOptions) -> NamedTuple

Simulate Lagrangian particle transport for a snow crab larval cohort across nursery
benthic settlement filtering.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.

# Outputs
- `NamedTuple`: Trajectories record and summary metrics.
"""
function run_segment_tracking(; opts::HydrodynamicOptions = HydrodynamicOptions(),
                              hydro_file::Union{Nothing, AbstractString} = nothing)
    println("\n=================================================================")
    println(" [Segment 6/8] Lagrangian Particle Tracking & Larval Behavior")
    println("=================================================================")
    rng = MersenneTwister(opts.seed)

    # The realized ocean simulation is a Tethys store, and the coupling between the two
    # halves is that store's identity -- a path. TethysFlowRun.jl prints it as `--hydro-file`,
    # which the CLI passes here on `opts.hydro_model_file`; nothing else carries it across, and
    # the ocean half's own configuration is never opened here. Its interior cell centres also
    # state the domain, so the spawn area and the advective field agree by construction.
    hydro_span = [0.0, 0.0]
    # When overwrite is requested, purge existing trajectory checkpoints and Zarr store
    if opts.overwrite
        tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
        for ext in (".zarr", ".jld2")
            tagged_cp = joinpath(opts.output_dir, "larval_trajectories$(tag)$(ext)")
            default_cp = joinpath(opts.output_dir, "larval_trajectories$(ext)")
            ispath(tagged_cp) && rm(tagged_cp; recursive = true, force = true)
            ispath(default_cp) && rm(default_cp; recursive = true, force = true)
        end
        println("Flag --overwrite active: purged existing trajectory checkpoints in $(opts.output_dir)")
    end

    # The store identity: the `--hydro-file` path when passed, else the run's default store.
    location = (hydro_file !== nothing && !isempty(strip(String(hydro_file)))) ?
        abspath(String(hydro_file)) : resolve_hydro_model_path(opts, "hydrodynamics.zarr")[1]
    store = try
        open_hydrodynamics(location)
    catch err
        if opts.track_only
            error("--track-only needs a realized simulation to track through; none is at " *
                  "$(location) ($(sprint(showerror, err))). Run TethysFlowRun.jl --hydro-only " *
                  "first, or pass --hydro-file=<store>.")
        end
        @warn "No Tethys store at $(location) ($(sprint(showerror, err))); " *
              "the analytical shelf jet will stand in."
        nothing
    end

    # Spawning area: the CFA boundaries if present, else the realized simulation's own domain
    # (its interior cell centres, so larvae are released where the flow is defined), else the
    # configured domain expanded by the buffer. Intersect with the hydrodynamic store bounds
    # when a store is loaded so larvae are never seeded outside the realized ocean grid.
    cfa_polys = load_cfa_polygons(opts.input_dir)
    spawn_lon, spawn_lat = if !isempty(cfa_polys)
        env = get_strata_buffered_envelope(cfa_polys, buffer_km = opts.buffer_km)
        (env.lon_range, env.lat_range)
    elseif !isnothing(store)
        (extrema(store.axes.lon), extrema(store.axes.lat))
    else
        Tethys.expand_domain_with_buffer(opts.domain_lon, opts.domain_lat, buffer_km = opts.buffer_km)
    end
    if !isnothing(store)
        s_lon = extrema(store.axes.lon)
        s_lat = extrema(store.axes.lat)
        spawn_lon = (max(spawn_lon[1], s_lon[1]), min(spawn_lon[2], s_lon[2]))
        spawn_lat = (max(spawn_lat[1], s_lat[1]), min(spawn_lat[2], s_lat[2]))
    end

    # Resolve coastline polygons BEFORE spawning so land rejection uses real shoreline
    # geometry at both initialization and transport steps.
    coast_candidates = String[
        joinpath(opts.input_dir, "coastline.parquet"),
        isnothing(store) ? "" : joinpath(dirname(store.location), "inputs", "coastline.parquet"),
        joinpath(Tethys.store_path("natural_earth_10m"), "coastline.parquet")
    ]
    coast_path = findfirst(p -> !isempty(p) && isfile(p), coast_candidates)
    coast_polys = isnothing(coast_path) ? nothing :
                  Tethys.load_coastline_polygons(coast_candidates[coast_path])
    if !isnothing(coast_polys) && isempty(coast_polys)
        coast_polys = nothing
    end
    if !isnothing(coast_polys)
        println("Coastline: loaded $(length(coast_polys)) polygon(s) from $(coast_candidates[coast_path]).")
    end

    # Seabed bathymetric elevation. Built *before* anything uses it, so the spawn and
    # the advection sit on the same sea floor.
    run_bathy = _larval_bathymetry(opts, store)
    bathy_field_fn = get_bathymetry_interpolator(
        run_bathy.elevation, Float64.(collect(run_bathy.lon)), Float64.(collect(run_bathy.lat)))

    # Build the hydrodynamic flow interpolator BEFORE spawning so `wet_mask_fn` can reject
    # cells that lie outside the realized (possibly rotated) ocean footprint.
    opts.hydro_start_days < 0.0 && error(
        "hydro_start_days must be non-negative (got $(opts.hydro_start_days))."
    )
    flow_interpolator = if !isnothing(store)
        start_time_req = store.times[1] + opts.hydro_start_days * 86400.0
        println("Reading 4D simulated fields from the Tethys store at $(store.location) " *
                "($(length(store.times)) interval(s), start offset = $(opts.hydro_start_days) d)...")
        try
            interp = build_flow_from_store(
                store;
                max_snapshots = opts.max_flow_snapshots,
                start_time = start_time_req,
                span = hydro_span
            )
            println("  retained span: $(round(hydro_span[1], digits = 1)) .. " *
                    "$(round(hydro_span[2], digits = 1)) s " *
                    "($(round((hydro_span[2] - hydro_span[1]) / 86400, digits = 2)) days)")
            interp
        catch err
            if opts.allow_analytical_fallback
                @warn "Failed to read the flow field from $(location): " *
                      "$(sprint(showerror, err)). Falling back to the analytical jet because " *
                      "allow_analytical_fallback is set -- results will NOT reflect the simulation."
                nothing
            else
                error(
                    "Failed to read the flow field from $(location):\n  " *
                    sprint(showerror, err) * "\n\n" *
                    "Tracking needs the simulated fields; substituting an analytical current would " *
                    "silently change the science. To proceed deliberately with a synthetic flow, " *
                    "re-run in analytic-fallback mode. To reduce memory, lower " *
                    "[hydrodynamics].max_flow_snapshots (currently $(opts.max_flow_snapshots))."
                )
            end
        end
    else
        nothing
    end

    # Align tracking clock t = 0 with the retained hydrodynamic origin hydro_span[1]
    # (which includes any TethysFlow spinup offset plus opts.hydro_start_days).
    # Build a static 3D BitArray wet mask ONCE from the initial hydrodynamic slice +
    # coastline polygons so wet/land checks take O(log N) nanoseconds and never evict
    # the active two-snapshot sliding window.
    hydro_t0 = hydro_span[1]
    domain_interior_fn = if !isnothing(flow_interpolator) && !isnothing(store)
        _build_domain_interior_fn(flow_interpolator; margin_km = 25.0, rot_deg = 35.0)
    else
        nothing
    end
    wet_mask_fn = if !isnothing(flow_interpolator) && !isnothing(store)
        raw_wet = _build_static_wet_mask(
            flow_interpolator;
            bathy_elevation = run_bathy.elevation,
            coastline = coast_polys
        )
        raw_wet
    else
        nothing
    end
    spawn_wet_fn = if !isnothing(wet_mask_fn) && !isnothing(domain_interior_fn)
        (lo, la, z) -> wet_mask_fn(lo, la, z) && domain_interior_fn(lo, la)
    else
        wet_mask_fn
    end

    println("Initializing $(opts.n_particles) Zoea I larvae in marine water (mode = $(opts.release_depth_mode), min depth >= $(opts.min_seabed_depth) m)...")
    larvae = initialize_larval_particles(
        opts.n_particles,
        lon_range = spawn_lon,
        lat_range = spawn_lat,
        release_depth_mode = opts.release_depth_mode,
        bottom_offset = opts.bottom_release_offset,
        ascent_target_depth = opts.ascent_target_depth,
        min_seabed_depth = opts.min_seabed_depth,
        max_seabed_depth = opts.max_seabed_depth,
        buffer_km = 0.0, # Buffer is already incorporated into spawn_lon/spawn_lat
        bathymetry = bathy_field_fn,
        coastline = coast_polys,
        wet_mask_fn = spawn_wet_fn,
        stage = :zoea1,
        rng = rng
    )

    flow_field_fn = if !isnothing(flow_interpolator)
        fw = flow_interpolator
        t0 = hydro_t0
        (lon, lat, z, t) -> begin
            f = fw(lon, lat, z, t + t0)
            (f.u, f.v, f.w)
        end
    else
        # Scotian Shelf alongshore current: flows southwestward along the shelf edge.
        # u ≈ -0.08 to -0.12 m/s (westward), v ≈ -0.03 to -0.06 m/s (southward/offshore).
        # Reference: Loder, J. W. & Petrie, B. (1991), CJFAS.
        (lon, lat, z, t) -> begin
            depth_decay  = exp(z / 120.0)          # velocity decays with depth
            lon_norm = clamp((lon - opts.domain_lon[1]) /
                             (opts.domain_lon[2] - opts.domain_lon[1]), 0.0, 1.0)
            jet_factor = 0.5 + 0.8 * exp(-((lon_norm - 0.3)^2) / 0.04)
            # The mean flow only. There is no tidal term here: tide is carried by the
            # ocean half's stored flow field, and an analytical fallback has no store to
            # carry it. A sinusoid added here would be a tide the model never ran.
            u_mean = -0.10 * jet_factor * depth_decay
            v_mean = -0.04 * jet_factor * depth_decay
            w_mean = 0.00020 * sin(2.0 * π * lon_norm) * depth_decay
            (u_mean, v_mean, w_mean)
        end
    end

    temp_field_fn = if !isnothing(flow_interpolator)
        fw = flow_interpolator
        t0 = hydro_t0
        (lon, lat, z, t) -> _sample_temperature(fw, lon, lat, z, t + t0)
    else
        # 3D Scotian Shelf thermal structure: surface mixed layer + CIL + slope water.
        # Consistent with set_initial_stratification! in hydrodynamic_model.jl.
        (lon, lat, z, t) -> begin
            lon_min, lon_max = opts.domain_lon
            lat_min, lat_max = opts.domain_lat
            x_norm = clamp((lon - lon_min) / (lon_max - lon_min), 0.0, 1.0)
            y_norm = clamp((lat - lat_min) / (lat_max - lat_min), 0.0, 1.0)
            T_surf   = 14.0 + 8.0 * x_norm - 5.0 * y_norm   # surface mixed layer
            T_cil_min = 1.5
            T_slope   = 8.5
            if z > -20.0
                frac = (z + 20.0) / 20.0
                clamp(T_cil_min + frac * (T_surf - T_cil_min), -1.5, 22.0)
            elseif z > -80.0
                centre_frac = (z + 50.0) / 30.0
                clamp(T_cil_min + 2.5 * centre_frac^2, -1.5, T_surf)
            else
                depth_factor = (abs(z) - 80.0) / 100.0
                clamp(T_cil_min + depth_factor * (T_slope - T_cil_min), T_cil_min, T_slope)
            end
        end
    end

    # Honor the configured track_duration_days (<= 0 uses the full remaining hydrodynamic
    # span), capped to the available hydrodynamic span when a hydrodynamic store is loaded.
    # The cohort still stops early once every larva is dead or settled.
    track_horizon = if !isnothing(flow_interpolator) && hydro_span[2] > hydro_span[1]
        avail_span = hydro_span[2] - hydro_span[1]
        if opts.track_duration <= 0.0
            avail_span
        else
            if opts.track_duration > avail_span
                @warn "Configured track_duration_days ($(opts.track_duration / 86400.0)) exceeds " *
                      "available hydrodynamic span ($(avail_span / 86400.0) days); capping to available span."
            end
            min(opts.track_duration, avail_span)
        end
    else
        opts.track_duration
    end
    # Probabilistic egg hatching aligned with TethysBio's temperature-dependent brooding
    # rule (`tau = 365 * (5 / max(1, T_bed))^0.3 + N(0, 5)`), integrated against the
    # realized local benthic temperature climatology when the tracking horizon covers
    # at least half an annual cycle.
    release_times = if track_horizon >= 180.0 * 86400.0
        base_date = Dates.Date(2020, 6, 1)  # Default June 1, 2020 origin matching scotian_shelf.toml
        if !isnothing(store) && hasproperty(store, :attrs)
            iso_str = get(store.attrs, "time_origin_iso", "")
            if !isempty(string(iso_str))
                try
                    base_date = Dates.Date(first(split(string(iso_str), "T")))
                catch
                end
            end
        end
        nom_start_date = base_date + Dates.Day(round(Int, hydro_t0 / 86400.0))
        origin_doy = Float64(Dates.dayofyear(nom_start_date))
        rt = sample_brood_hatch_times(
            larvae.lon, larvae.lat, larvae.depth;
            temperature_fn = temp_field_fn,
            origin_doy = origin_doy,
            total_duration = track_horizon,
            rng = rng
        )
        println("  probabilistic hatching (Arrhenius-weighted Jan-Apr cohort, start $(Dates.format(nom_start_date, "yyyy-mm-dd")), DOY $(round(Int, origin_doy))): " *
                "release days $(round(minimum(rt) / 86400.0 + 1.0, digits=1))..$(round(maximum(rt) / 86400.0 + 1.0, digits=1)) d " *
                "(median $(round(median(rt) / 86400.0 + 1.0, digits=1)) d)")
        rt
    else
        nothing
    end

    println("Tracking cohort over up to $(track_horizon / 86400.0) days " *
            "(dt=$(opts.track_dt)s, stops early when all larvae are dead or settled)...")
    trajectories = track_larval_cohort(
        larvae,
        velocity_fn = flow_field_fn,
        temperature_fn = temp_field_fn,
        bathymetry_fn = bathy_field_fn,
        total_duration = track_horizon,
        dt = opts.track_dt,
        output_interval = opts.output_schedule_seconds,
        release_times = release_times,
        max_current_speed = opts.max_current_speed,
        κ_h = opts.diffusivity_h,
        κ_v = opts.diffusivity_v,
        is_lat_lon = true,
        enable_molting = opts.enable_molting,
        # Stochasticity controls. These must be passed explicitly: `track_larval_cohort` defaults
        # all three coefficients of variation to 0.0, i.e. a fully deterministic cohort. Loading
        # them into `opts` and writing them to `resolved_config.toml` is not enough on its own --
        # if they are not forwarded here the provenance file records a stochastic model that the
        # run never actually used.
        cv_molt = opts.cv_molt,
        cv_mortality = opts.cv_mortality,
        cv_settlement = opts.cv_settlement,
        settlement_stochastic = opts.settlement_stochastic,
        # The biology's parameters, forwarded from the configuration. The behaviour module
        # has its own defaults for all of these; forwarding them is what makes the
        # configuration file the record of the model that actually ran.
        temp_base_degreedays = opts.temp_base_degreedays,
        degree_days_larval = opts.degree_days_larval,
        intermolt_days_larval = opts.intermolt_days_larval,
        larval_baseline_hazards = opts.larval_baseline_hazards,
        larval_temp_opt = opts.larval_temp_opt,
        larval_temp_sigma = opts.larval_temp_sigma,
        larval_overage_hazard = opts.larval_overage_hazard,
        larval_max_intermolt_days = opts.larval_max_intermolt_days,
        settlement_min_depth = opts.settlement_min_depth,
        settlement_max_depth = opts.settlement_max_depth,
        settlement_max_temp = opts.settlement_max_temp,
        zoea1_day_depth = opts.zoea1_day_depth,
        zoea1_night_depth = opts.zoea1_night_depth,
        zoea2_day_depth = opts.zoea2_day_depth,
        zoea2_night_depth = opts.zoea2_night_depth,
        megalopa_day_depth = opts.megalopa_day_depth,
        megalopa_night_depth = opts.megalopa_night_depth,
        zoea2_depth_factor = opts.zoea2_depth_factor,
        megalopa_swim_factor = opts.megalopa_swim_factor,
        enable_bbl = true,
        enable_sinking = true,
        enable_initial_ascent = opts.enable_initial_ascent,
        ascent_speed = opts.ascent_speed,
        ascent_target_depth = opts.ascent_target_depth,
        # When `wet_mask_fn` is active it already encodes the baked coastline polygons;
        # avoid repeating polygon ray-casting on every particle step.
        coastline = isnothing(wet_mask_fn) ? coast_polys : nothing,
        wet_mask_fn = wet_mask_fn,
        domain_interior_fn = domain_interior_fn,
        rng = rng
    )

    # Release the streaming hydrodynamic slice cache before writing the Zarr checkpoint.
    flow_interpolator = nothing
    flow_field_fn = nothing
    temp_field_fn = nothing
    wet_mask_fn = nothing

    # Save trajectories to Zarr checkpoint for modular reloading
    mkpath(opts.output_dir)
    tag = !isempty(opts.run_id) ? "_$(opts.run_id)" : ""
    track_checkpoint = joinpath(opts.output_dir, "larval_trajectories$(tag).zarr")
    _to_f32(x) = isnothing(x) ? nothing : (x isa AbstractArray{Float64} ? Float32.(x) : x)
    cp_payload = Dict{String, Any}(
        "lons" => _to_f32(trajectories.lons),
        "lats" => _to_f32(trajectories.lats),
        "depths" => _to_f32(trajectories.depths),
        "temperatures" => hasproperty(trajectories, :temperatures) ? _to_f32(trajectories.temperatures) : nothing,
        "degree_days" => _to_f32(trajectories.degree_days),
        "degree_days_timeseries" => hasproperty(trajectories, :degree_days_timeseries) ? _to_f32(trajectories.degree_days_timeseries) : nothing,
        "survival_probability" => hasproperty(trajectories, :survival_probability) ? _to_f32(trajectories.survival_probability) : nothing,
        "stage_survival" => hasproperty(trajectories, :stage_survival) ? _to_f32(trajectories.stage_survival) : nothing,
        "stages" => trajectories.stages,
        "alive" => trajectories.alive,
        "settlement_status" => trajectories.settlement_status,
        "settlement_age" => _to_f32(trajectories.settlement_age),
        "ascent_duration" => hasproperty(trajectories, :ascent_duration) ? _to_f32(trajectories.ascent_duration) : nothing,
        # The cohort-level stage diagnostics are saved too, so `--viz` can redraw the
        # progression figures from a checkpoint without re-running Segment 6. `molt_schedule`
        # is a tuple of `MoltCDF` structs and is reconstructed deterministically from `opts`
        # by `_complete_stage_diagnostics`.
        "cohort_molt_fraction" => hasproperty(trajectories, :cohort_molt_fraction) ? _to_f32(trajectories.cohort_molt_fraction) : nothing,
        "molt_schedule" => nothing,
        "times" => trajectories.times,
        "ids" => trajectories.ids,
    )
    tmp_checkpoint = track_checkpoint * ".tmp"
    ispath(tmp_checkpoint) && rm(tmp_checkpoint; recursive = true, force = true)
    Tethys.write_data(tmp_checkpoint, cp_payload; backend = :zarr, overwrite = true)
    ispath(track_checkpoint) && rm(track_checkpoint; recursive = true, force = true)
    mv(tmp_checkpoint, track_checkpoint; force = true)

    # Maintain default checkpoint as fallback when custom run_id is supplied
    if !isempty(strip(opts.run_id))
        default_cp = joinpath(opts.output_dir, "larval_trajectories.zarr")
        try
            isdir(default_cp) && rm(default_cp; recursive = true, force = true)
            cp(track_checkpoint, default_cp; force = true)
        catch
        end
    end

    # Save trajectories directly into Zarr
    if opts.enable_zarr
        try
            println("Archiving trajectories to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
            try
                target_run_id = run_identifier(opts)
                save_simulation_run!(
                    db, target_run_id, opts;
                    trajectories = trajectories,
                    config = options_to_configuration(opts),
                    notes = "Larval tracking run" *
                            (!isempty(strip(opts.run_id)) ? " [$(opts.run_id)]" : "")
                )
                println("Trajectories for '$(target_run_id)' successfully archived in Zarr.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive trajectories to Zarr: $(err)"
        end
    end

    n_settled = count(==( :settled_successful), trajectories.settlement_status)
    n_alive   = count(identity, trajectories.alive)
    println("Lagrangian tracking complete:")
    println("  Total particles:     $(opts.n_particles)")
    println("  Surviving particles: $(n_alive) / $(opts.n_particles)")
    println("  Settled on nursery:  $(n_settled) / $(opts.n_particles) " *
            "($(round(100.0 * n_settled / opts.n_particles, digits=1))%)")
    println("  Checkpoint saved:    $(track_checkpoint)")

    return (trajectories = trajectories, checkpoint_path = track_checkpoint)
end


"""
    run_segment_metrics(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Compute gridded larval retention and recruitment metrics, thermal exposure indices,
empirical Lagrangian drift velocities and effective diffusivities, regional
connectivity matrices across Crab Fishing Areas (CFAs), and export multi-layer NetCDF and Zarr.

# Inputs
- `opts::HydrodynamicOptions`: Workflow options.
- `trajectories`: Optional trajectory NamedTuple.

# Outputs
- `NamedTuple`: Metric outputs, connectivity matrix, and export filepaths.
"""
function run_segment_metrics(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    # The tracking segment returns a *record of a run* -- the cohort plus where it was
    # checkpointed -- while every consumer below works on the cohort itself. Unwrapping once
    # here, rather than in each consumer, keeps the two records distinct. Passing the run
    # record is what produced Trajectory structure must contain :lons and :lats from
    # inside the analysis, which named the wrong layer.
    raw_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories
    target_trajs = hasproperty(raw_trajs, :trajectories) ? raw_trajs.trajectories : raw_trajs

    println("\n=================================================================")
    println(" [Segment 7/8] Empirical Movement, Recruitment & Connectivity")
    println("=================================================================")

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    # 1. Empirical Advection & Turbulent Diffusivity
    println("Estimating empirical advection and turbulent diffusivity fields...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 2. Gridded Recruitment & Settlement Metrics
    println("Computing gridded recruitment and benthic nursery settlement metrics...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 3. Gridded Thermal Exposure & Degree-Days
    println("Computing gridded thermal exposure metrics...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)

    # 4. Demographic Connectivity Matrix across Crab Fishing Areas (CFAs)
    println("Computing macro-regional demographic connectivity matrix...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    cfa_definitions = if !isempty(cfa_polys)
        println("Loaded $(length(cfa_polys)) CFA boundary polygons from $(opts.input_dir)/: $(join([p.name for p in cfa_polys], ", "))")
        vcat(cfa_polys, [(name = "Outside / Offshore",)])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_definitions)
    println("Transition Probability Matrix (P_ij):")
    for i in 1:length(conn.strata_names)
        row_str = join([string(round(conn.matrix[i, j], digits = 3)) for j in 1:length(conn.strata_names)], ", ")
        println("  $(rpad(conn.strata_names[i], 28)) -> [$(row_str)]")
    end

    # 5. Export multi-variable NetCDF4 and Zarr archives
    nc_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.nc")
    jld_export_path = joinpath(opts.output_dir, "larval_dispersal_analysis.zarr")
    active_config = options_to_configuration(opts)

    println("Exporting multi-layer NetCDF4 archive -> $(nc_export_path)...")
    export_larval_dispersal_netcdf(
        nc_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

    println("Exporting Zarr archive -> $(jld_export_path)...")
    export_larval_dispersal_jld2(
        jld_export_path,
        trajectories = target_trajs,
        lon_bins = lon_b,
        lat_bins = lat_b,
        strata_definitions = cfa_definitions,
        config = active_config
    )

# 7. Archive simulation run, trajectories, metrics & connectivity in Zarr
    if opts.enable_zarr
        try
            println("Archiving simulation run and metrics to Zarr -> $(opts.zarr_path)...")
            db = open_storage(opts.zarr_path)
            try
                target_run_id = run_identifier(opts)
                save_simulation_run!(
                    db,
                    target_run_id,
                    opts;
                    trajectories = target_trajs,
                    metrics = (
                        mean_exposure_temperature = let vals = filter(!isnan, vec(therm_metrics.mean_exposure_temperature))
                            isempty(vals) ? 0.0 : Statistics.mean(vals)
                        end,
                        mean_degree_days = let vals = filter(!isnan, vec(therm_metrics.mean_degree_days))
                            isempty(vals) ? 0.0 : Statistics.mean(vals)
                        end,
                        settlement_success_rate = count(==(:settled_successful), target_trajs.settlement_status) /
                            max(1, length(target_trajs.settlement_status))
                    ),
                    connectivity = conn,
                    gridded_dispersal = (
                        lon_centers = emp_mov.lon_centers,
                        lat_centers = emp_mov.lat_centers,
                        u_mean = emp_mov.u_mean,
                        v_mean = emp_mov.v_mean,
                        diffusivity = emp_mov.diffusivity,
                        density = rec_metrics.settlement_density,
                        mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                        mean_degree_days = therm_metrics.mean_degree_days,
                        sample_count = emp_mov.sample_count
                    ),
                    config = active_config,
                    notes = "Larval workflow run" *
                            (!isempty(strip(opts.run_id)) ? " [$(opts.run_id)]" : "")
                )
                println("Zarr run '$(target_run_id)' successfully archived.")
            finally
                close_storage(db)
            end
        catch err
            @warn "Failed to archive simulation run and metrics to Zarr: $(err)"
        end
    end

    return (
        empirical_movement = emp_mov,
        recruitment_metrics = rec_metrics,
        thermal_metrics = therm_metrics,
        connectivity = conn,
        netcdf_path = nc_export_path,
        jld2_path = jld_export_path,
        storage_path = opts.enable_zarr ? opts.zarr_path : nothing
    )
end


"""
    render_all_figures(; opts) -> NamedTuple

Regenerate every figure and visual for an existing run, from stored results only.

This is the implementation behind `--figures`: particle figures, the larval-biology diagnostics,
the Eulerian hydrodynamic set, the interactive map, and the animation. It reads the trajectory
checkpoint and the hydrodynamic archive rather than re-simulating, so it is the fast path to
iterate on a plotting function.

# What it draws
- **Larval** (7): trajectories, DVM profiles, settlement density, empirical movement field,
  connectivity matrix, thermal exposure, recruitment summary, plus the three biology diagnostics
  (molt progression, degree-day growth, survival curves).
- **Eulerian** (9): advection, tracers, stratification, diffusion, section, station time series,
  time–depth diagram, near-bed temperature map, T–S diagram.
- **Interactive**: the decimated HTML dashboard.
- **Animation**: the field/dashboard animation, when `--animate-hydro` is set.

Missing inputs are reported rather than fatal, so a run with, say, no hydrodynamic archive still
produces the larval figures.
"""
function render_all_figures(; opts::HydrodynamicOptions = HydrodynamicOptions())
    println("\n=================================================================")
    println(" [Figures] Regenerating all figures and visuals")
    println("=================================================================")
    return run_segment_visualize(opts = opts)
end


"""
    _ck_or(saved::Dict, key::AbstractString, default)

Read `key` from a loaded JLD2 checkpoint dictionary, falling back to `default` when the key is
absent **or** stored as `nothing`.

The tracking checkpoint writes optional fields as `nothing` when the trajectory lacks them, and
older checkpoints predate several fields entirely, so a plain `get(saved, key, default)` is not
enough — a stored `nothing` would flow into the figures as a missing value.
"""
_ck_or(saved, key::AbstractString, default) =
    (haskey(saved, key) && !isnothing(saved[key])) ? saved[key] : default


"""
    _trajectories_from_checkpoint(saved::Dict) -> NamedTuple

Build the trajectory NamedTuple from a loaded `larval_trajectories*.zarr` dictionary.

Fields that older checkpoints lack are filled with the same placeholders used before
(`temperatures` 4 °C, `degree_days_timeseries` 40 °C·d, `survival_probability` 0.95). The
stage diagnostics `cohort_molt_fraction` and `molt_schedule` are set to `nothing` when absent
and are then rebuilt by [`_complete_stage_diagnostics`](@ref); a zero placeholder would plot
as a cohort that never moults.
"""
function _trajectories_from_checkpoint(saved)
    n_p, n_t = size(saved["lons"])
    t_end = saved["times"][end]
    return (
        lons = saved["lons"],
        lats = saved["lats"],
        depths = saved["depths"],
        temperatures = _ck_or(saved, "temperatures", fill(4.0, n_p, n_t)),
        degree_days = saved["degree_days"],
        degree_days_timeseries = _ck_or(saved, "degree_days_timeseries",
                                        fill(40.0, n_p, n_t)),
        survival_probability = _ck_or(saved, "survival_probability", fill(0.95, n_p, n_t)),
        stages = Symbol.(saved["stages"]),
        alive = Bool.(saved["alive"]),
        settlement_status = Symbol.(saved["settlement_status"]),
        settlement_age = _ck_or(saved, "settlement_age", fill(t_end, n_p)),
        stage_survival = _ck_or(saved, "stage_survival", fill(1.0, n_p)),
        ascent_duration = _ck_or(saved, "ascent_duration", fill(0.0, n_p)),
        cohort_molt_fraction = _ck_or(saved, "cohort_molt_fraction", nothing),
        molt_schedule = _ck_or(saved, "molt_schedule", nothing),
        times = saved["times"],
        ids = saved["ids"]
    )
end


"""
    _complete_stage_diagnostics(trajs::NamedTuple, opts) -> NamedTuple

Add `molt_schedule` and `cohort_molt_fraction` to a trajectory set that lacks them, such as one
read back from Zarr/Tethys storage, which stores per-particle rows only.

Both are deterministic functions of the run configuration and the stored trajectory, so they
are recomputed rather than approximated:

- `molt_schedule = molt_schedule(degree_days_larval = opts.degree_days_larval,
  cv_molt = opts.cv_molt)`. The schedule is rebuilt from the CONFIGURED thresholds, so
  the diagnostic agrees with the run it is diagnosing. It used to be rebuilt with the
  package's defaults, which made the figure a picture of a different model whenever the
  configuration said anything else -- and the configuration is where the model's
  parameters live, so it nearly always says something else.
- `cohort_molt_fraction[:, j]` is the fraction of the cohort that has completed each
  transition at step `j`, computed from the stored stage history. It used to be computed
  from the cohort-MEAN degree-day, which was meaningful when the thresholds were
  cumulative and every larva's accumulation counted from the same origin. They are now
  accumulated PER STAGE and reset at each ecdysis, so the cohort mean is an average over
  larvae in different stages -- each counting from its own last ecdysis -- and reading
  the threshold CDF at that mean says more about the stage mix than about the cohort's
  development. With the three per-stage thresholds equal, it also saturates the second
  and third transitions while every larva is still in the first.

The stage history is the source of truth for this quantity, and it carries the dead
larvae's stages forward: a larva recorded as `:dead` is counted at the stage it reached.

Fields already present are left unchanged. When molting is disabled the run recorded no
schedule, so nothing is added and the figure reports it as not recorded.
"""
function _complete_stage_diagnostics(trajs::NamedTuple, opts)
    has(k) = hasproperty(trajs, k) && !isnothing(getproperty(trajs, k))
    opts.enable_molting || return trajs

    sched = has(:molt_schedule) ? trajs.molt_schedule :
             molt_schedule(degree_days_larval = opts.degree_days_larval, cv_molt = opts.cv_molt)
    extra = has(:molt_schedule) ? NamedTuple() : (molt_schedule = sched,)

    if !has(:cohort_molt_fraction)
        stages = getproperty(trajs, :stages)
        extra = merge(extra, (cohort_molt_fraction = cohort_transition_fraction_history(stages),))
    end
    return isempty(extra) ? trajs : merge(trajs, extra)
end


"""
    load_run_trajectories(opts::HydrodynamicOptions) -> NamedTuple

Load the trajectory set of the configured run for Segments 7 and 8.

Sources, in order:
1. Analytical storage (`opts.zarr_path`, run id `run_identifier(opts)`). If the Zarr
   checkpoint describes the same cohort (equal `ids` and `times`), the checkpoint is used
   instead, because it also holds `stage_survival`, `ascent_duration` and the true
   `settlement_age`, none of which analytical storage stores.
2. The Zarr checkpoint `larval_trajectories[_<run_id>].zarr`.
3. Neither present: Segment 6 is run.

The result is passed through [`_complete_stage_diagnostics`](@ref), so every source yields the
same set of stage fields.
"""
function load_run_trajectories(opts::HydrodynamicOptions)
    target_run_id = run_identifier(opts)
    tag = !isempty(strip(opts.run_id)) ? "_$(opts.run_id)" : ""
    tagged_cp = joinpath(opts.output_dir, "larval_trajectories$(tag).zarr")
    track_cp = isdir(tagged_cp) ? tagged_cp :
               joinpath(opts.output_dir, "larval_trajectories.zarr")

    from_db = nothing
    if !opts.overwrite && opts.enable_zarr && isdir(opts.zarr_path)
        db = open_storage(opts.zarr_path; read_only = true)
        try
            println("Loading particle trajectories from Zarr for run '$(target_run_id)'...")
            from_db = load_trajectories_namedtuple(db, target_run_id)
        catch err
            println("  Zarr load failed ($(sprint(showerror, err))); trying the checkpoint.")
        finally
            close_storage(db)
        end
    end

    from_cp = if isdir(track_cp)
        try
            raw_cp = Tethys.read_data(track_cp; backend = :zarr)
            if haskey(raw_cp, "lons") && haskey(raw_cp, "lats") && haskey(raw_cp, "times")
                _trajectories_from_checkpoint(raw_cp)
            else
                @warn "Checkpoint at $(track_cp) is incomplete (missing 'lons'/'lats'/'times'); removing and re-running Segment 6."
                rm(track_cp; recursive = true, force = true)
                nothing
            end
        catch err
            @warn "Could not read checkpoint at $(track_cp) ($(sprint(showerror, err))); removing and re-running Segment 6."
            rm(track_cp; recursive = true, force = true)
            nothing
        end
    else
        nothing
    end

    trajs = if !isnothing(from_db) && !isnothing(from_cp) &&
               Int.(from_cp.ids) == Int.(from_db.ids) &&
               Float64.(from_cp.times) == Float64.(from_db.times)
        println("  using the matching checkpoint $(track_cp) (holds the full stage record).")
        from_cp
    elseif !isnothing(from_db)
        from_db
    elseif !isnothing(from_cp)
        println("Loading particle trajectories from checkpoint: $(track_cp)...")
        from_cp
    else
        println("Trajectories checkpoint not found. Running Segment 6...")
        run_segment_tracking(opts = opts).trajectories
    end
    return _complete_stage_diagnostics(trajs, opts)
end


"""
    run_segment_visualize(;
        opts::HydrodynamicOptions,
        trajectories=nothing
    ) -> NamedTuple

Generate and export publication-ready CairoMakie figures and standalone
    interactive Leaflet HTML dashboards with multi-run comparison and
rich spatial layers.

# Outputs
- `NamedTuple`: Paths to generated PNG figures and HTML dashboard.
"""
function run_segment_visualize(;
    opts::HydrodynamicOptions = HydrodynamicOptions(),
    trajectories = nothing
)
    # The tracking segment returns a *record of a run* -- the cohort plus where it was
    # checkpointed -- while every consumer below works on the cohort itself. Unwrapping once
    # here, rather than in each consumer, keeps the two records distinct. Passing the run
    # record is what produced Trajectory structure must contain :lons and :lats from
    # inside the analysis, which named the wrong layer.
    raw_trajs = isnothing(trajectories) ? load_run_trajectories(opts) : trajectories
    target_trajs = hasproperty(raw_trajs, :trajectories) ? raw_trajs.trajectories : raw_trajs

    println("\n=================================================================")
    println(" [Segment 8/8] Scientific Visualizations & Spatial Figures")
    println("=================================================================")
    mkpath(opts.output_dir)

 

    # Bathymetry for the background contours. Open the configured hydrodynamic store if
    # present so `_larval_bathymetry` reads the store's own `bottom_height` relief first.
    hydro_location = resolve_hydro_model_path(opts, "hydrodynamics.zarr")[1]
    hydro_store = isdir(hydro_location) ?
        (try open_hydrodynamics(hydro_location) catch; nothing end) : nothing
    bathy_data = _larval_bathymetry(opts, hydro_store)

    viz_coast_candidates = String[
        joinpath(opts.input_dir, "coastline.parquet"),
        isnothing(hydro_store) ? "" : joinpath(dirname(hydro_store.location), "inputs", "coastline.parquet"),
        joinpath(Tethys.store_path("natural_earth_10m"), "coastline.parquet")
    ]
    viz_coast_idx = findfirst(p -> !isempty(p) && isfile(p), viz_coast_candidates)
    viz_coast = isnothing(viz_coast_idx) ? nothing :
                Tethys.load_coastline_polygons(viz_coast_candidates[viz_coast_idx])
    if !isnothing(viz_coast) && isempty(viz_coast)
        viz_coast = nothing
    end

    lon_b = range(opts.domain_lon[1], opts.domain_lon[2], length = opts.grid_size[1])
    lat_b = range(opts.domain_lat[1], opts.domain_lat[2], length = opts.grid_size[2])

    map_rot_deg = 35.0

    # 1. 2D Particle Trajectories Map
    fig1_path = joinpath(opts.output_dir, "larval_trajectories.png")
    println("Rendering 2D spatial trajectories map -> $(fig1_path)...")
    cfa_polys = load_cfa_polygons(opts.input_dir)
    plot_particle_trajectories(
        target_trajs,
        bathymetry_data = bathy_data,
        strata = cfa_polys,
        coastline = viz_coast,
        rotation_deg = map_rot_deg,
        title = "Snow Crab Larval Dispersal",
        output_path = fig1_path
    )

    # 2. DVM Depth Profiles
    fig2_path = joinpath(opts.output_dir, "dvm_depth_profiles.png")
    println("Rendering DVM depth profiles -> $(fig2_path)...")
    plot_vertical_migration_profiles(
        target_trajs,
        sample_indices = 1:min(5, size(target_trajs.depths, 1)),
        title = "Larval Diel Vertical Migration (DVM) Profiles",
        output_path = fig2_path
    )

    # 3. 2D Settlement Nursery Density Heatmap
    fig3_path = joinpath(opts.output_dir, "settlement_density.png")
    println("Rendering 2D settlement density distribution -> $(fig3_path)...")
    plot_larval_dispersal_density(
        target_trajs,
        coastline = viz_coast,
        rotation_deg = map_rot_deg,
        title = "Snow Crab Nursery Settlement Density (%)",
        output_path = fig3_path
    )

    # 4. Empirical Movement Velocity Field and Diffusivity
    fig4_path = joinpath(opts.output_dir, "empirical_movement_field.png")
    println("Rendering empirical velocity quiver field -> $(fig4_path)...")
    emp_mov = estimate_empirical_movement(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_empirical_movement_field(
        emp_mov;
        coastline = viz_coast,
        rotation_deg = map_rot_deg,
        output_path = fig4_path
    )

    # 5. Annotated Regional Connectivity Matrix
    fig5_path = joinpath(opts.output_dir, "regional_connectivity_matrix.png")
    println("Rendering regional connectivity matrix -> $(fig5_path)...")
    cfa_defs = if !isempty(cfa_polys)
        vcat(cfa_polys, [(name = "Outside / Offshore",)])
    else
        [
            (name = "CFA 20-22 (Eastern NS)", lon = (-62.0, -57.0), lat = (44.5, 47.5)),
            (name = "CFA 23-24 (Middle Shelf)", lon = (-64.5, -60.0), lat = (43.0, 45.5)),
            (name = "CFA 4X (Southwest NS)", lon = (-68.0, -64.0), lat = (42.0, 44.5)),
            (name = "Offshore / Slope", lon = (-68.0, -57.0), lat = (40.0, 43.0))
        ]
    end
    conn = compute_empirical_connectivity(target_trajs, strata_definitions = cfa_defs)
    plot_connectivity_matrix(conn, output_path = fig5_path)

    # 6. Thermal Exposure Map
    fig6_path = joinpath(opts.output_dir, "thermal_exposure_map.png")
    println("Rendering thermal exposure map -> $(fig6_path)...")
    therm_metrics = compute_gridded_thermal_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_thermal_exposure_map(
        therm_metrics;
        coastline = viz_coast,
        rotation_deg = map_rot_deg,
        output_path = fig6_path
    )

    # 7. Recruitment Summary
    fig7_path = joinpath(opts.output_dir, "recruitment_summary.png")
    println("Rendering recruitment summary bar chart -> $(fig7_path)...")
    rec_metrics = compute_gridded_recruitment_metrics(target_trajs, lon_bins = lon_b, lat_bins = lat_b)
    plot_recruitment_summary(
        rec_metrics;
        coastline = viz_coast,
        rotation_deg = map_rot_deg,
        output_path = fig7_path
    )

    # 7b. Larval biology diagnostics: stage progression, thermal-time growth, and survival.
    # These were previously only embedded in the interactive HTML payload and never drawn, which
    # left the quantile/CDF stage model unverifiable from the run output.
    fig_molt_path = joinpath(opts.output_dir, "molt_progression.png")
    fig_dd_path = joinpath(opts.output_dir, "degree_day_growth.png")
    fig_surv_path = joinpath(opts.output_dir, "survival_curves.png")
    println("Rendering molt progression -> $(fig_molt_path)...")
    plot_molt_progression(target_trajs, output_path = fig_molt_path)
    println("Rendering degree-day growth -> $(fig_dd_path)...")
    plot_degree_day_growth(target_trajs, output_path = fig_dd_path)
    println("Rendering survival curves -> $(fig_surv_path)...")
    plot_survival_curves(target_trajs, output_path = fig_surv_path)

    # 8. Hydrodynamic Model Eulerian Advection & Tracers Figures
    fig8_path = joinpath(opts.output_dir, "hydrodynamic_advection.png")
    fig9_path = joinpath(opts.output_dir, "hydrodynamic_tracers.png")
    fig_strat_path = joinpath(opts.output_dir, "hydrodynamic_stratification.png")
    fig_diff_path = joinpath(opts.output_dir, "hydrodynamic_diffusion.png")
    fig_sec_path = joinpath(opts.output_dir, "hydrodynamic_section.png")
    fig_ts_path = nothing
    fig_hov_path = nothing
    fig_tbot_path = nothing
    fig_ts_diagram_path = nothing

    # The Eulerian figures are drawn from the REALIZED OCEAN STORE -- the same store the
    # cohort tracked through -- and not from the larval Zarr store's `hydrodynamic_fields`
    # group. That group is created empty by the schema and never written, so reading it
    # returned an empty dataset: no depth axis, no speed, no time, and the first figure
    # died on `depths[0]`. The ocean store is where the Eulerian fields actually are, and
    # `extract_hydrodynamic_dataset` returns them in the shape these figures take.
    hydro_location = resolve_hydro_model_path(opts, "hydrodynamics.zarr")[1]
    active_hydro = if isdir(hydro_location)
        try
            extract_hydrodynamic_dataset(hydro_location)
        catch err
            @warn "Could not read the hydrodynamic fields from $(hydro_location) " *
                  "($(sprint(showerror, err))); the Eulerian figures are skipped."
            nothing
        end
    else
        @warn "No realized ocean store at $(hydro_location); the Eulerian figures are skipped."
        nothing
    end

    if !isnothing(active_hydro)
        println("Rendering hydrodynamic advection velocity field -> $(fig8_path)...")
        plot_hydrodynamic_advection(
            active_hydro;
            bathymetry_data = bathy_data,
            coastline = viz_coast,
            rotation_deg = map_rot_deg,
            output_path = fig8_path
        )

        println("Rendering hydrodynamic seawater temperature & salinity tracers -> $(fig9_path)...")
        plot_hydrodynamic_tracers(
            active_hydro;
            coastline = viz_coast,
            rotation_deg = map_rot_deg,
            output_path = fig9_path
        )

        println("Rendering hydrodynamic stratification diagnostics (N², S) -> $(fig_strat_path)...")
        plot_hydrodynamic_stratification(
            active_hydro;
            coastline = viz_coast,
            rotation_deg = map_rot_deg,
            output_path = fig_strat_path
        )

        # Skip vertical diffusivity map when the store carries only a uniform background placeholder.
        diff_vals = filter(isfinite, vec(active_hydro.diffusion))
        if !isempty(diff_vals) && (maximum(diff_vals) - minimum(diff_vals)) > 1e-12
            println("Rendering hydrodynamic turbulent diffusion & eddy viscosity -> $(fig_diff_path)...")
            plot_hydrodynamic_diffusion(
                active_hydro;
                coastline = viz_coast,
                rotation_deg = map_rot_deg,
                output_path = fig_diff_path
            )
        else
            println("Skipping vertical diffusion figure (κ_v is not archived in the hydrodynamic store).")
        end

        println("Rendering across-shelf hydrodynamic vertical cross-section -> $(fig_sec_path)...")
        plot_hydrodynamic_section(
            active_hydro;
            variable = :temperature,
            section_type = :cross_shelf,
            output_path = fig_sec_path
        )

        # 8b. Eulerian time series at a representative wet shelf station.
        fig_ts_path = joinpath(opts.output_dir, "hydrodynamic_timeseries.png")
        println("Rendering hydrodynamic station time series -> $(fig_ts_path)...")
        try
            plot_hydrodynamic_timeseries(
                active_hydro;
                store = hydro_store,
                station = (-62.0, 44.0),
                variable = :temperature,
                output_path = fig_ts_path
            )
        catch err
            @warn "Station time-series plot failed (non-fatal): $(sprint(showerror, err))"
        end

        # 8c. Additional Eulerian diagnostics.
        fig_hov_path = joinpath(opts.output_dir, "hydrodynamic_depth_time_profile.png")
        fig_tbot_path = joinpath(opts.output_dir, "bottom_temperature_map.png")
        fig_ts_diagram_path = joinpath(opts.output_dir, "ts_diagram.png")
        st_lon, st_lat = -62.0, 44.0
        for (label, path, f) in (
            ("depth-time profile", fig_hov_path,
             () -> plot_hydrodynamic_hovmoller(active_hydro; store = hydro_store,
                                               station = (st_lon, st_lat),
                                               variable = :temperature, output_path = fig_hov_path)),
            ("near-bed temperature map", fig_tbot_path,
             () -> plot_hydrodynamic_bottom_temperature(active_hydro; coastline = viz_coast,
                                                        rotation_deg = map_rot_deg,
                                                        output_path = fig_tbot_path)),
            ("temperature-salinity diagram", fig_ts_diagram_path,
             () -> plot_hydrodynamic_ts_diagram(active_hydro; output_path = fig_ts_diagram_path))
            )
            println("Rendering $label -> $(path)...")
            try
                f()
            catch err
                # A diagnostic that cannot be built from this archive must not abort the run.
                @warn "$label failed (non-fatal): $(first(split(sprint(showerror, err), "\n")))"
            end
        end
    end


    # 9. Hydrodynamic Field & Dashboard Animation (Oceananigans / CairoMakie)
    anim_file_path = nothing
    if opts.animate_hydro
        anim_fmt = lowercase(opts.anim_format) == "gif" ? "gif" : "mp4"
        hydro_source = if !isnothing(active_hydro)
            active_hydro
        else
            default_store = "hydrodynamics.zarr"
            resolved_store, _ = resolve_hydro_model_path(opts, default_store)
            ispath(resolved_store) ? resolved_store : nothing
        end

        anim_file_path = if !isempty(strip(opts.anim_output_path))
            opts.anim_output_path
        elseif opts.anim_variable == :dashboard
            joinpath(opts.output_dir, "hydrodynamic_dashboard_animation.$(anim_fmt)")
        else
            joinpath(opts.output_dir, "hydrodynamic_$(opts.anim_variable)_animation.$(anim_fmt)")
        end

        if opts.anim_variable == :dashboard
            println("Rendering animated hydrodynamic dashboard -> $(anim_file_path)...")
            animate_hydrodynamic_dashboard(
                hydro_source;
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        else
            println("Rendering animated hydrodynamic $(opts.anim_variable) field -> $(anim_file_path)...")
            animate_hydrodynamic_field(
                hydro_source;
                variable = opts.anim_variable,
                depth = opts.anim_depth,
                trajectories = target_trajs,
                bathymetry_data = bathy_data,
                output_path = anim_file_path,
                framerate = opts.anim_fps,
                show_trajectories = opts.anim_overlay_particles,
                domain_lon = opts.domain_lon,
                domain_lat = opts.domain_lat
            )
        end
    end

    # 9. Query Archived Runs from Zarr for Cross-Run Comparison & Multi-Layer Interactive Map
    runs_bundle = Dict{String, Any}()
    if opts.enable_zarr && isdir(opts.zarr_path)
        try
            db = open_storage(opts.zarr_path; read_only = true)
            try
                runs_df = list_simulation_runs(db)
                for r in eachrow(runs_df)
                    r_id = string(r.run_id)
                    try
                        r_trajs = load_trajectories_namedtuple(db, r_id)
                        r_disp = try load_gridded_dispersal(db, r_id) catch; nothing end
                        r_conn = try load_connectivity_matrix(db, r_id) catch; nothing end
                        r_hydro = try load_hydrodynamic_field(db, r_id) catch; nothing end
                        runs_bundle[r_id] = (
                            trajectories = r_trajs,
                            gridded_dispersal = r_disp,
                            connectivity = r_conn,
                            hydrodynamics = r_hydro
                        )
                    catch
                    end
                end
            finally
                close_storage(db)
            end
        catch
        end
    end

    if isempty(runs_bundle) && !isnothing(target_trajs)
        runs_bundle[run_identifier(opts)] = (
            trajectories = target_trajs,
            gridded_dispersal = (
                lon_centers = emp_mov.lon_centers,
                lat_centers = emp_mov.lat_centers,
                u_mean = emp_mov.u_mean,
                v_mean = emp_mov.v_mean,
                diffusivity = emp_mov.diffusivity,
                density = rec_metrics.settlement_density,
                mean_exposure_temperature = therm_metrics.mean_exposure_temperature,
                mean_degree_days = therm_metrics.mean_degree_days,
                sample_count = emp_mov.sample_count
            ),
            connectivity = conn,
            hydrodynamics = active_hydro
        )
    end

    # Cross-run 2D figure
    fig10_path = joinpath(opts.output_dir, "run_comparison.png")
    println("Rendering cross-run comparison -> $(fig10_path)...")
    run_comp = Dict{Symbol, Any}()
    for (r_label, r_data) in runs_bundle
        run_comp[Symbol(r_label)] = r_data.trajectories
    end

    compare_run_dispersal(
        run_comp,
        title = "Snow Crab Larval Dispersal Across Archived Runs",
        output_path = fig10_path
    )

    # 10. Interactive HTML5 Dashboard (with zlib-compressed sidecar payload)
    html_path = joinpath(opts.output_dir, "interactive_larval_tracks.html")
    if opts.interactive_map
        base_date = Dates.Date(2020, 6, 1)
        if !isnothing(hydro_store) && hasproperty(hydro_store, :attrs)
            iso_str = get(hydro_store.attrs, "time_origin_iso", "")
            if !isempty(string(iso_str))
                try
                    base_date = Dates.Date(first(split(string(iso_str), "T")))
                catch
                end
            end
        end
        store_t0_days = (!isnothing(hydro_store) && !isempty(hydro_store.times)) ?
            (Float64(hydro_store.times[1]) / 86400.0) : 0.0
        nom_start_date = base_date + Dates.Day(round(Int, store_t0_days + opts.hydro_start_days))
        start_date_iso = Dates.format(nom_start_date, "yyyy-mm-dd")
        println("Rendering interactive Leaflet larval dispersal dashboard -> $(html_path) (start date: $(start_date_iso))...")
        export_interactive_tracks_html(
            html_path;
            runs_data = runs_bundle,
            title = "Scotian Shelf Snow Crab Larval Dispersal Trajectories",
            start_date_iso = start_date_iso,
            max_points_per_track = 550
        )
    end

    println("All visualization figures and interactive maps successfully generated.")
    return (
        fig_trajectories = fig1_path,
        fig_dvm = fig2_path,
        fig_density = fig3_path,
        fig_empirical_movement = fig4_path,
        fig_connectivity = fig5_path,
        fig_thermal = fig6_path,
        fig_recruitment = fig7_path,
        # Larval biology diagnostics
        fig_molt_progression = fig_molt_path,
        fig_degree_day_growth = fig_dd_path,
        fig_survival = fig_surv_path,
        # Eulerian
        fig_hydro_advection = fig8_path,
        fig_hydro_tracers = fig9_path,
        fig_hydro_stratification = fig_strat_path,
        fig_hydro_diffusion = fig_diff_path,
        fig_hydro_section = fig_sec_path,
        fig_hydro_timeseries = fig_ts_path,
        fig_hovmoller = fig_hov_path,
        fig_bottom_temperature = fig_tbot_path,
        fig_ts_diagram = fig_ts_diagram_path,
        fig_comparison = fig10_path,
        interactive_map = opts.interactive_map ? html_path : nothing,
        animation = anim_file_path
    )
end


"""
    run_cli_list_runs(; opts::HydrodynamicOptions)

Query and display all simulation runs currently archived in Zarr.
"""
function run_cli_list_runs(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Archived Simulation Runs in Zarr: $(opts.zarr_path)")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = list_simulation_runs(db)
    close_storage(db)

    if nrow(df) == 0
        println("No simulation runs archived in database yet.")
        return
    end

    for row in eachrow(df)
        println("-----------------------------------------------------------------")
        println("Run ID:          $(row.run_id)")
        println("Created:         $(row.created_at)")
        println("Particles:       $(row.n_particles) larvae | Duration: $(round(row.duration_days, digits=1)) days")
        println("Recruitment:     $(round(row.settlement_success_rate * 100, digits=2))% successful")
        println("Mean PLD:        $(round(row.mean_pld_days, digits=1)) days | Mean Temp: $(round(row.mean_exposure_temperature, digits=2)) °C")
        println("Dispersal Dist:  $(round(row.mean_dispersal_distance_km, digits=1)) km")
    end
    println("-----------------------------------------------------------------\n")
    return df
end


"""
    run_cli_compare_runs(; opts::HydrodynamicOptions)

Query and print comparative metrics across archived runs in Zarr.
"""
function run_cli_compare_runs(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Run Comparative Analysis from Zarr")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    df = compare_runs(db)
    close_storage(db)

    if nrow(df) == 0
        println("No simulation data available for run comparison.")
        return
    end

    for row in eachrow(df)
        println("Run: $(rpad(row.run_id, 16)) | Archived runs: $(row.n_runs)")
        println("  Settlement Success: $(round(row.mean_settlement_success * 100, digits=2))%")
        println("  Mean PLD:           $(round(row.mean_pld_days, digits=1)) days")
        println("  Mean Exposure Temp: $(round(row.mean_temperature_celsius, digits=2)) °C")
        println("  Mean Dispersal:     $(round(row.mean_dispersal_km, digits=1)) km")
        println("  Thermal Mortality:  $(round(row.mean_thermal_mortality_rate * 100, digits=2))%")
        println()
    end
    return df
end


"""
    run_cli_model_average(; opts::HydrodynamicOptions)

Compute and display ensemble model-averaged demographic connectivity and recruitment.
"""
function run_cli_model_average(; opts::HydrodynamicOptions = HydrodynamicOptions())
    if !isdir(opts.zarr_path)
        println("No Zarr storage found at: $(opts.zarr_path)")
        return
    end
    println("\n=================================================================")
    println(" Multi-Run Ensemble Model Averaging")
    println("=================================================================")
    db = open_storage(opts.zarr_path; read_only = true)
    runs_df = list_simulation_runs(db)

    if nrow(runs_df) == 0
        close_storage(db)
        println("No simulation data available for model averaging.")
        return
    end

    run_ids = unique(runs_df.run_id)
    println("Averaging across $(length(run_ids)) archived runs: $(join(run_ids, ", "))")
    ens = compute_ensemble_model_average(db, run_ids)
    close_storage(db)

    println("\nEnsemble Weighted Settlement Success: $(round(ens.mean_recruitment_rate * 100, digits=2))%")
    println("Ensemble Weighted Mean PLD:          $(round(ens.mean_pld_days, digits=1)) days")
    println("Ensemble Weighted Thermal Exposure:  $(round(ens.mean_thermal_exposure, digits=2)) °C")
    println("\nEnsemble Model-Averaged Connectivity Matrix (P_ij ± std):")
    n_s = length(ens.strata_names)
    for i in 1:n_s
        row_strs = [
            "$(round(ens.mean_connectivity[i, j], digits=3))±$(round(ens.std_connectivity[i, j], digits=3))"
            for j in 1:n_s
        ]
        println("  $(rpad(ens.strata_names[i], 28)) -> [$(join(row_strs, ", "))]")
    end
    println()
    return ens
end

"""
    animate_hydrodynamic_dashboard(args...; kwargs...)

Owned by the OCEAN half. TethysDrift calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
animate_hydrodynamic_dashboard(args...; kwargs...) = error("`animate_hydrodynamic_dashboard` is owned by TethysFlow, not by this package.")

"""
    animate_hydrodynamic_field(args...; kwargs...)

Owned by the OCEAN half. TethysDrift calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
animate_hydrodynamic_field(args...; kwargs...) = error("`animate_hydrodynamic_field` is owned by TethysFlow, not by this package.")

"""
    _larval_bathymetry(opts, store)

The seabed `(lon, lat) -> z` the larvae sit on, resolved without calling the ocean half.

Order, and why each step is what it is:

1. **The store's own `bottom_height`**, from its `.zattrs`. The run that wrote the store put the
   seabed it immersed *into* the store, so a larva released where the store says it is land is
   measured against the sea floor the flow was computed on. This is the coupling the two halves
   agreed on, and it is a path, not a package import.
2. **The data layer**, for the store's own domain. `Tethys` is shared by both halves and
   is not the ocean half: it resolves a bathymetry product by coverage out of the central
   lakehouse. A store written before the seabed was carried in `.zattrs` has no step 1, and this
   is the alternative that answers with the same product the ocean run would have resolved.

Nothing is invented. A store with no seabed and a lakehouse with no product is an error, rather
than a flat or analytic sea floor that would spawn larvae in the air or on a ridge the run never
used.
"""
function _larval_bathymetry(opts::HydrodynamicOptions, store)
    # 1. The store's own seabed.
    if store !== nothing
        zattrs = joinpath(store.location, ".zattrs")
        if isfile(zattrs)
            a = JSON3.read(read(zattrs, String))
            if haskey(a, "bottom_height") && haskey(store.axes, :lon)
                lon = collect(Float64, store.axes.lon)
                lat = collect(Float64, store.axes.lat)
                bh = _parse_bottom_height(a, length(lon), length(lat))
                size(bh) == (length(lon), length(lat)) || error(
                    "The store at $(store.location) carries a bottom_height of size " *
                    "$(size(bh)) for a $(length(lon)) x $(length(lat)) grid. The seabed a " *
                    "run immersed must be the shape of that grid; refusing the read.")
                println("Seabed: the store's own immersed relief, $(size(bh)).")
                return (elevation = bh, lon = lon, lat = lat)
            end
        end
    end

    # 2. The data layer, for the domain the store states (or the configured one without a store).
    if store !== nothing
        lon_r = (minimum(store.axes.lon), maximum(store.axes.lon))
        lat_r = (minimum(store.axes.lat), maximum(store.axes.lat))
    else
        lon_r, lat_r = opts.domain_lon, opts.domain_lat
    end
    println("Seabed: no relief in the store; resolving a bathymetry product from the central " *
            "lakehouse, as the ocean half would.")
    # A durable destination in the lakehouse, keyed by the resolved source. Without one the
    # data layer writes to `inputs/bathymetry.zarr` in the caller's working directory -- which
    # is this package's directory, is not shared with the ocean half, and cannot be rewritten
    # once it exists. Naming the store makes the product machine-shared, as every other
    # assimilated input is, and makes a second run a cache hit.
    src = Tethys.preferred_bathymetry_source(; lon_range = lon_r, lat_range = lat_r)
    destination = joinpath(Tethys.store_path(String(src)), "bathymetry.zarr")
    arr = try
        ds = if ispath(destination)
            Tethys.load_dataset(destination)
        else
            Tethys.bathymetry_dataset(;
                source = src, lon_range = lon_r, lat_range = lat_r,
                output_path = destination, verbose = false)
        end
        elev_var = haskey(ds.variables, "elevation") ? "elevation" :
                   first(Tethys.variables_like(ds, ["elevation", "altitude", "z", "topo"]))
        elev_ga = ds.variables[elev_var]
        perm = Tethys.dim_permutation(elev_ga.dims, :lon, :lat)
        lon_c, _ = Tethys.axis(ds, :lon)
        lat_c, _ = Tethys.axis(ds, :lat)
        raw = Array{Float64}(elev_ga.data)
        elev = ndims(raw) == 1 ? reshape(raw, :, 1) : permutedims(raw, perm)
        (elevation = elev, lon = lon_c, lat = lat_c)
    catch err
        error("No seabed is available. The store carries no `bottom_height` and the central " *
              "lakehouse holds no bathymetry product for lon $(lon_r) lat $(lat_r) " *
              "($(sprint(showerror, err))). A larval run placed on a sea floor nobody " *
              "measured would be a different simulation: run TethysFlowRun.jl once to " *
              "assimilate a product, or fetch one with `Tethys.fetch_input`.")
    end
    println("  seabed: $(length(arr.lon)) x $(length(arr.lat)) grid, source $(src)")
    return (elevation = Float64.(arr.elevation),
            lon = Float64.(collect(arr.lon)),
            lat = Float64.(collect(arr.lat)))
end

"""
    _run_bathymetry_source(args...; kwargs...)

Owned by the OCEAN half. TethysDrift calls into neither package's ocean half: importing it would be the edge SPLIT_PLAN.md section 1 forbids. A larval `--viz` that renders a hydrodynamic field is this package describing the ocean, and should read the store the ocean half wrote.
"""
_run_bathymetry_source(args...; kwargs...) = error("`_run_bathymetry_source` is owned by TethysFlow, not by this package.")
