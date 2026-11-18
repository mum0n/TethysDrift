# ==============================================================================
#  store_reader.jl -- reading the realized simulation from its Tethys store
# ==============================================================================
#
#  Replaces `flow_reader.jl`, which read Oceananigans' JLD2 archive plus a separate grid
#  sidecar. That reader existed because the archive was the only thing a run produced; now the
#  run writes its slices into a Tethys store that carries its own coordinates, so a reader
#  opens one path and gets flow, axes, and times from it.
#
#  What the store provides (written by `TethysFlow.attach_output_callback!`):
#
#      hydrodynamics.zarr/          a TethysDataset
#        .zattrs                    title, created_utc, and the seabed the run immersed
#        lon, lat, depth, time      coordinate arrays, in the dataset's own `coords`
#        T, u, v, w                 fields, each (lon, lat, depth, time)
#
#  Reading it is `Tethys.load_dataset(store)`. That is the whole point of the store being a Tethys
#  store: an earlier writer put each interval in its own directory (`T/0`, `T/1`), which Tethys's
#  Zarr reader cannot open -- it discovers variables as the root group's child arrays -- so every
#  consumer walked Zarr itself. This did. Now it does not.
#
#  This package does not know how the store was produced. It knows the layout, and that is the
#  whole coupling between the two packages: `TethysFlow` writes it, this package reads it,
#  and neither imports the other. The layout is asserted in `test_store_reader.jl` rather than
#  assumed, so a change on either side shows as a failed read.
#
#  What is deliberately NOT here: checkpoints. Restart is the ocean half's concern and stays
#  in JLD2, the one format Oceananigans can checkpoint-restart from.
# ==============================================================================

using Tethys
using Tethys: TethysDataset

"""
    _store_path(location) -> String

Resolve a location to the store directory.

A directory is the store itself; a path whose suffix is `.zarr` is used as given. Anything else
is an error naming what was received rather than a guess at what was meant, because a wrong
guess here reads a different simulation than the one asked for.
"""
function _store_path(location::AbstractString)
    p = abspath(location)
    isdir(p) && return p
    endswith(p, ".zarr") && isdir(p) && return p
    error("No hydrodynamic store at $(p). Expected a directory, or a `.zarr` store. This is " *
          "the location the run was handed, and a missing store cannot be substituted.")
end

"""
    _read_times(store) -> Vector{Float64}

The time of each written interval, in order, from the store's own `time` coordinate.

One axis for the whole store is deliberate: a per-field copy could drift, and a field read at
its own recorded time but indexed by a shared time axis is a silent mismatch. The coordinate
is read here from a dataset opened with `load_dataset`, so it is the same axis the ocean half wrote.
"""
function _read_times(store::AbstractString)
    ds = _load_hydrodynamic_dataset(store)
    haskey(ds.coords, :time) || error(
        "The store at $(store) has no `time` coordinate, so the times of its intervals are " *
        "unknown. A field with no time is not a simulation output.")
    out = collect(Float64, ds.coords[:time].data)
    isempty(out) && error("The store at $(store) records no times.")
    return out
end

"""
    _read_axes(store) -> (; lon, lat, depth)

The grid's interior cell centres, from the store's own dimension records.

Present for the reason the sidecar used to be: these are what a reader needs to place a value
on Earth, and a store that does not state cannot be interpolated against. They come from the
dataset rather than from a second file that nothing kept in step, and a store that declares
fewer than three axes errors naming which, rather than defaulting to a spanned range --
defaulting would put every value in a track at the wrong place.
"""
function _read_axes(store::AbstractString)
    ds = _load_hydrodynamic_dataset(store)
    out = (lon = nothing, lat = nothing, depth = nothing)
    for (k, name) in ((:lon, "longitude"), (:lat, "latitude"), (:depth, "depth"))
        haskey(ds.dims, k) || haskey(ds.dims, Symbol(name)) || error(
            "The store at $(store) declares $(sort!(collect(keys(ds.dims)))) but not $(k). " *
            "A store without its cell centres cannot be read.")
        c = haskey(ds.dims, k) ? collect(Float64, ds.dims[k].coords) :
                                 collect(Float64, ds.dims[Symbol(name)].coords)
        out = merge(out, NamedTuple{(k,)}((c,)))
    end
    return out
end

"""
    hydrodynamic_variables(store) -> Vector{String}

Which fields the store carries. A reader that asked for a field the store does not have would
otherwise return zeros and call it a velocity.
"""
function hydrodynamic_variables(store::AbstractString)
    return sort!(collect(keys(_load_hydrodynamic_dataset(store).variables)))
end

"""
    _load_hydrodynamic_dataset(store) -> TethysDataset

The store as a Tethys dataset, opened through Tethys.

An empty variable set is an error here rather than at the first read: it is what a store that
is not Tethys's shape produces -- `load_dataset` succeeds and finds nothing -- and the failure
used to surface three call sites later as "does not carry T", naming the wrong layer.
"""
function _load_hydrodynamic_dataset(store::AbstractString)
    ds = Tethys.load_dataset(_store_path(store); backend = :zarr)
    isempty(ds.variables) && error(
        "The store at $(store) opens as a Tethys dataset with no variables. Either it is " *
        "empty, or it was written by an older writer that put each interval in its own " *
        "directory (`T/0`, `T/1`), which Tethys cannot read. Re-run the simulation to " *
        "produce a Tethys store.")
    return ds
end

"""
    open_hydrodynamics(location; fields = nothing) -> HydrodynamicStore

Open a realized simulation's store for reading.

`fields` defaults to everything the store carries, so a caller gets the whole run rather than
only what it happened to name. A request for a field the store does not have is an error
listing what it does have.

The result holds the opened dataset. The intervals are read on demand, which is what lets a
consumer hold a long run without holding it in memory.
"""
function open_hydrodynamics(location::AbstractString; fields = nothing)
    store = _store_path(location)
    ds = _load_hydrodynamic_dataset(store)
    times = _read_times(store)
    axes = _read_axes(store)
    available = sort!(collect(keys(ds.variables)))

    wanted = isnothing(fields) ? available :
             [String(f) for f in (fields isa Union{Tuple, AbstractVector} ? fields : (fields,))]
    unknown = [f for f in wanted if !(f in available)]
    isempty(unknown) || error(
        "The store at $(store) does not carry $(unknown). It carries $(available).")

    return HydrodynamicStore(store, available, times, axes, ds)
end

"""
    HydrodynamicStore

A realized simulation's store, opened for reading: its location, its fields, its times, its
axes, and the dataset they came from. No flow data.

The separation is the point. A consumer is handed a location and reads the interval it needs,
rather than being handed every interval the run produced.
"""
struct HydrodynamicStore
    location::String
    variables::Vector{String}
    times::Vector{Float64}
    axes::@NamedTuple{lon::Vector{Float64}, lat::Vector{Float64},
                      depth::Vector{Float64}}
    dataset::TethysDataset
end

"""
    variable(store, name, m) -> Array{Float64, 3}

Interval `m` of one field, as `(n_lon, n_lat, n_depth)`.

Read from the dataset's own array, at `[:, :, :, m]`. An out-of-range index is an error naming
the interval count, not a clamped read: a track that asked for interval 7 of a 4-interval run
is a bug in the track, and returning the last interval would hide it behind a plausible
velocity.
"""
function variable(store::HydrodynamicStore, name::AbstractString, m::Integer)
    !(name in store.variables) && error(
        "The store at $(store.location) does not carry `$(name)`. It carries " *
        "$(store.variables).")
    (1 <= m <= length(store.times)) || error(
        "Interval $(m) is out of range: this run wrote $(length(store.times)) interval(s).")
    v = store.dataset.variables[name]
    size(v.data, 4) >= m || error(
        "Interval $(m) is declared by the time axis but field `$(name)` has only " *
        "$(size(v.data, 4)) interval(s). The store and its time axis disagree.")
    # Copied to the host rather than returned as a disk-backed view: the interpolator reads
    # this slice many times per tracking step, and a view is a decompression each time.
    return Array{Float64}(v.data[:, :, :, m])
end

"""
    create_flow_interpolator(location; fields = ("u", "v", "w", "T"))

Build the `(lon, lat, z, t) -> (u, v, w, T)` function the larval model integrates against.

This replaces `create_flow_interpolator_from_jld2`. The difference is where the axes come from:
they are read from the store rather than re-derived from a sidecar, so a reader cannot disagree
with the run about where its cells are.

Bilinear in space and linear in time, with the nearest edge value outside the domain. A track
that wanders past the boundary still returns a real velocity rather than failing, which is
what `track_larval_cohort` expects of a flow field.
"""
function create_flow_interpolator(location::AbstractString;
                                  fields::Tuple = ("u", "v", "w", "T"))
    store = open_hydrodynamics(location; fields = fields)
    lon = store.axes.lon
    lat = store.axes.lat
    z = store.axes.depth
    t = store.times

    nx, ny, nz = length(lon), length(lat), length(z)
    nt = length(t)

    # Hold the arrays for the run. The store is read here rather than per step, because
    # `track_larval_cohort` samples it many times per interval, and re-reading a Zarr array
    # per sample turns the tracking step into I/O.
    data = Dict{String, Vector{Array{Float64, 3}}}()
    for name in fields
        data[name] = [variable(store, name, m) for m in 1:nt]
    end

    # The trilinear weights are the same for every field, so they are computed once per call
    # rather than once per field per step.
    _bracket(vals, v) = begin
        n = length(vals)
        n <= 1 && return (1, 1, 0.0)
        i = clamp(something(findfirst(x -> x > v, vals), n), 2, n)
        denom = vals[i] - vals[i - 1]
        s = denom == 0.0 ? 0.0 : clamp((v - vals[i - 1]) / denom, 0.0, 1.0)
        return (i - 1, i, s)
    end

    return function flow_interpolator(lon_q::Real, lat_q::Real, z_q::Real, t_q::Real)
        i0, i1, sx = _bracket(lon, Float64(lon_q))
        j0, j1, sy = _bracket(lat, Float64(lat_q))
        k0, k1, sz = _bracket(z, Float64(z_q))
        m0, m1, θ = _bracket(t, Float64(t_q))

        function sample(name)
            arr = data[name]
            v0 = arr[m0]; v1 = arr[m1]
            a = 0.0
            for (ar, wf) in ((v0, 1.0 - θ), (v1, θ))
                a += wf * (
                    (1 - sx) * (1 - sy) * (1 - sz) * ar[i0, j0, k0] +
                    sx       * (1 - sy) * (1 - sz) * ar[i1, j0, k0] +
                    (1 - sx) * sy       * (1 - sz) * ar[i0, j1, k0] +
                    sx       * sy       * (1 - sz) * ar[i1, j1, k0] +
                    (1 - sx) * (1 - sy) * sz       * ar[i0, j0, k1] +
                    sx       * (1 - sy) * sz       * ar[i1, j0, k1] +
                    (1 - sx) * sy       * sz       * ar[i0, j1, k1] +
                    sx       * sy       * sz       * ar[i1, j1, k1])
            end
            return a
        end

        return (u = sample("u"), v = sample("v"), w = sample("w"), T = sample("T"))
    end
end

"""
    store_location(location::AbstractString) -> String

The location a consumer is handed. A path, and nothing else.
"""
store_location(location::AbstractString) = abspath(location)

"""
    _parse_bottom_height(attrs, nx::Int, ny::Int) -> Matrix{Float64}

Parse `bottom_height` from `.zattrs`, handling both 2D nested JSON arrays and
1D flattened column-major vectors paired with `bottom_height_shape`.
"""
function _parse_bottom_height(attrs, nx::Int, ny::Int)::Matrix{Float64}
    raw = attrs["bottom_height"]
    if !isempty(raw) && first(raw) isa AbstractVector
        return Float64.(reduce(hcat, raw))
    end
    vec_bh = Float64.(collect(raw))
    if haskey(attrs, "bottom_height_shape") && length(attrs["bottom_height_shape"]) == 2
        shp = (Int(attrs["bottom_height_shape"][1]), Int(attrs["bottom_height_shape"][2]))
        if length(vec_bh) == prod(shp)
            return reshape(vec_bh, shp)
        end
    end
    length(vec_bh) == nx * ny && return reshape(vec_bh, nx, ny)
    return reshape(vec_bh, 1, :)
end

"""
    extract_hydrodynamic_dataset(location::AbstractString; time_index::Int = 1, time_seconds::Union{Nothing, Real} = nothing,
                                  domain_lon::Tuple{<:Real, <:Real} = (-180.0, 180.0),
                                  domain_lat::Tuple{<:Real, <:Real} = (-90.0, 90.0),
                                  grid_size::Tuple{Int, Int} = (0, 0)) -> NamedTuple

Read one interval from a TethysFlow Zarr store and return a NamedTuple compatible with
the figure pipeline. Fields: lons, lats, depths, u, v, w, temperature, salinity,
bathymetry, speed, diffusion, stratification, time_seconds, depth_index, depth_m.

This is TethysDrift's local reader for the shared store layout. The two packages do not
import each other; each defines its own reader for the store that TethysFlow writes.

# Inputs
- `location`: Path to a Zarr store directory written by TethysFlow.
- `time_index`: 1-based snapshot index.
- `time_seconds`: If provided, the nearest interval to this time is selected (overrides
  `time_index`).
- `domain_lon`, `domain_lat`: Accepted for API compatibility with TethysFlow's reader;
  the Zarr store already carries its own coordinates, so these are not used for subsetting.
- `grid_size`: Accepted for API compatibility; the Zarr store carries its own resolution.

# Outputs
A NamedTuple with fields matching TethysFlow's `extract_hydrodynamic_dataset` output,
so figure families can consume either package's store reader interchangeably.
"""
extract_hydrodynamic_dataset(hydro::NamedTuple; kwargs...) = hydro

function extract_hydrodynamic_dataset(location::AbstractString;
                                       time_index::Union{Nothing, Int} = nothing,
                                       time_seconds::Union{Nothing, Real} = nothing,
                                       domain_lon::Tuple{<:Real, <:Real} = (-180.0, 180.0),
                                       domain_lat::Tuple{<:Real, <:Real} = (-90.0, 90.0),
                                       grid_size::Tuple{Int, Int} = (0, 0))
    store = open_hydrodynamics(location)

    n_t = length(store.times)

    # Resolve time_index from time_seconds if provided
    if !isnothing(time_seconds)
        times = Float64.(store.times)
        time_index = argmin(abs.(times .- Float64(time_seconds)))
    elseif isnothing(time_index)
        time_index = 1
    end
    !(1 <= time_index <= n_t) && error(
        "Interval $time_index is out of range: this run wrote $(n_t) interval(s).")

    lons = store.axes.lon
    lats = store.axes.lat
    depths = store.axes.depth
    t_sec = Float64(store.times[time_index])

    # Read all fields at this interval
    u = variable(store, "u", time_index)
    v = variable(store, "v", time_index)
    w = variable(store, "w", time_index)
    T = variable(store, "T", time_index)

    # Salinity may or may not be in the store; default to reference 35 PSU if absent
    available = Set(store.variables)
    S = if "S" in available || "salinity" in available
        name = "S" in available ? "S" : "salinity"
        variable(store, name, time_index)
    else
        fill(35.0, size(T))
    end

    speed = sqrt.(u.^2 .+ v.^2)

    # Vertical eddy diffusivity κ_v — the store carries u, v, w, T only; diffusivity is a
    # derived closure field. Use a background value consistent with the model's closure
    # (1e-4 m²/s is typical shelf-ocean background mixing).
    diffusion = fill(1.0e-4, size(T))

    # Buoyancy frequency squared N² from vertical temperature gradient (linear EOS):
    #   N² = -g/ρ₀ * ∂ρ/∂z = g * α * ∂T/∂z   (thermal expansion only)
    # depths are negative below surface; ∂T/∂z computed by centered differences.
    g = 9.81      # m/s²
    alpha = 1.7e-4 # K⁻¹ thermal expansion coefficient
    nz = length(depths)
    stratification = zeros(Float64, size(T))
    if nz >= 2
        for k in 1:nz
            if k == 1
                dT_dz = (T[:, :, 2] - T[:, :, 1]) / (depths[2] - depths[1])
            elseif k == nz
                dT_dz = (T[:, :, nz] - T[:, :, nz - 1]) / (depths[nz] - depths[nz - 1])
            else
                dT_dz = (T[:, :, k + 1] - T[:, :, k - 1]) / (depths[k + 1] - depths[k - 1])
            end
            stratification[:, :, k] = g * alpha .* dT_dz
        end
    end

    # Bathymetry from .zattrs bottom_height, or 2D store variable if available, else wet shelf
    zattrs_path = joinpath(store.location, ".zattrs")
    bathy = nothing
    if isfile(zattrs_path)
        attrs = JSON3.read(read(zattrs_path, String))
        if haskey(attrs, "bottom_height")
            bh = _parse_bottom_height(attrs, length(lons), length(lats))
            if size(bh) == (length(lons), length(lats))
                bathy = bh
            end
        end
    end
    if isnothing(bathy) && ("bathymetry" in available || "depth" in available)
        name = "bathymetry" in available ? "bathymetry" : "depth"
        v_bathy = store.dataset.variables[name]
        if ndims(v_bathy.data) == 2
            bathy = Array{Float64}(v_bathy.data)
        elseif ndims(v_bathy.data) >= 3
            bathy = Array{Float64}(v_bathy.data[:, :, 1])
        end
    end
    if isnothing(bathy)
        bathy = fill(-100.0, length(lons), length(lats))
    end

    return (
        lons = lons,
        lats = lats,
        depths = depths,
        u = u,
        v = v,
        w = w,
        temperature = T,
        salinity = S,
        speed = speed,
        diffusion = diffusion,
        stratification = stratification,
        bathymetry = bathy,
        elevation = bathy,   # alias for interactive HTML which reads .elevation
        time_seconds = t_sec,
        depth_index = 1,
        depth_m = Float64(depths[1]),
    )
end
