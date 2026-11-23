# ==============================================================================
# Analysis, Single-Layer Store, Checkpoint & Visualization Tests
# ==============================================================================

using Test
using TethysDrift
using Tethys
using Zarr

@testset "Analysis, Store Edge Cases & Visualization" begin
    @testset "Single-layer / single-snapshot Zarr store" begin
        mktempdir() do tmp
            zpath = joinpath(tmp, "single_layer.zarr")
            lons = [-64.0, -63.0, -62.0]
            lats = [44.0, 45.0]
            depths = [-25.0]
            times = [0.0]

            root = Zarr.zgroup(zpath)
            _zput(g, name, arr) = (za = Zarr.zcreate(eltype(arr), g, name, size(arr)...); za[:] = arr)
            _zput(root, "lon", lons)
            _zput(root, "lat", lats)
            _zput(root, "z", depths)
            _zput(root, "time", times)

            nx, ny, nz, nt = length(lons), length(lats), length(depths), length(times)
            u = fill(0.1, nx, ny, nz, nt)
            v = fill(-0.05, nx, ny, nz, nt)
            w = zeros(nx, ny, nz, nt)
            T = fill(4.5, nx, ny, nz, nt)
            S = fill(32.5, nx, ny, nz, nt)
            bathy = fill(-120.0, nx, ny)

            _zput(root, "u", u)
            _zput(root, "v", v)
            _zput(root, "w", w)
            _zput(root, "T", T)
            _zput(root, "S", S)
            _zput(root, "bottom_height", bathy)

            hydro = TethysDrift.open_hydrodynamics(zpath)
            @test length(hydro.times) == 1
            @test length(hydro.axes.depth) == 1

            interp = TethysDrift.create_flow_interpolator(zpath)
            res = interp(-63.0, 44.5, -25.0, 100.0)
            @test isapprox(res.u, 0.1; atol = 1e-6)
            @test isapprox(res.v, -0.05; atol = 1e-6)
            @test isapprox(res.T, 4.5; atol = 1e-6)

            ds = TethysDrift.extract_hydrodynamic_dataset(zpath)
            @test size(ds.temperature) == (nx, ny, nz)
            @test all(ds.stratification .== 0.0)
            @test TethysDrift.extract_hydrodynamic_dataset(ds) === ds

            # Verify flat .zattrs bottom_height + bottom_height_shape parsing
            attrs_flat = Dict{String, Any}(
                "bottom_height" => vec(bathy),
                "bottom_height_shape" => [nx, ny]
            )
            parsed_bh = TethysDrift._parse_bottom_height(attrs_flat, nx, ny)
            @test size(parsed_bh) == (nx, ny)
            @test parsed_bh == bathy
        end
    end

    @testset "Domain buffer, strata envelope & coastline" begin
        expanded = Tethys.expand_domain_with_buffer((-64.0, -62.0), (44.0, 46.0); buffer_km = 20.0)
        @test expanded[1][1] < -64.0 && expanded[1][2] > -62.0
        @test expanded[2][1] < 44.0 && expanded[2][2] > 46.0

        strata = [
            (name = "S1", lons = [-64.0, -63.0, -63.0, -64.0], lats = [44.0, 44.0, 45.0, 45.0])
        ]
        lon_b, lat_b = TethysDrift.get_strata_buffered_envelope(strata; buffer_km = 10.0)
        @test lon_b[1] < -64.0 && lon_b[2] > -63.0
        @test lat_b[1] < 44.0 && lat_b[2] > 45.0

        polys = Tethys.load_coastline_polygons()
        @test polys isa AbstractVector
    end

    @testset "Checkpoint deserialization converts Symbol & Bool" begin
        saved = Dict{String, Any}(
            "lons" => [-63.5 -63.4; -63.0 -62.9],
            "lats" => [44.5 44.6; 44.8 44.9],
            "depths" => [-10.0 -12.0; -20.0 -22.0],
            "degree_days" => [50.0 55.0; 120.0 125.0],
            "stages" => ["zoea_1", "megalopa"],
            "alive" => [true, false],
            "settlement_status" => ["active", "settled_successful"],
            "times" => [0.0, 3600.0],
            "ids" => [1, 2]
        )
        trajs = TethysDrift._trajectories_from_checkpoint(saved)
        @test trajs.stages == [:zoea_1, :megalopa]
        @test trajs.alive == Bool[true, false]
        @test trajs.settlement_status == [:active, :settled_successful]
    end

    @testset "Eulerian & larval visualization smoke tests" begin
        lons = [-64.0, -63.0, -62.0]
        lats = [44.0, 45.0, 46.0]
        depths = [-10.0, -50.0, -100.0]
        nx, ny, nz = length(lons), length(lats), length(depths)

        T = [4.0 + 0.1 * i - 0.02 * k for i in 1:nx, j in 1:ny, k in 1:nz]
        S = [32.0 + 0.2 * k for i in 1:nx, j in 1:ny, k in 1:nz]
        T[1, 1, :] .= NaN
        S[1, 1, :] .= NaN
        u = fill(0.15, nx, ny, nz)
        v = fill(-0.05, nx, ny, nz)
        w = zeros(nx, ny, nz)
        speed = hypot.(u, v)
        bathy = fill(-120.0, nx, ny)

        hydro_nt = (
            lons = lons,
            lats = lats,
            depths = depths,
            times = [0.0],
            temperature = T,
            salinity = S,
            u = u,
            v = v,
            w = w,
            speed = speed,
            stratification = fill(1e-5, nx, ny, nz),
            diffusion = fill(1e-4, nx, ny, nz),
            viscosity = fill(1e-4, nx, ny, nz),
            density = fill(1025.0, nx, ny, nz),
            vorticity = zeros(nx, ny, nz),
            bathymetry = bathy,
            time_seconds = 0.0,
            n_snapshots = 1
        )

        @test TethysDrift.plot_hydrodynamic_tracers(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_stratification(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_diffusion(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_bottom_temperature_map(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_temperature_salinity_diagram(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_hovmoller(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_section(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_bottom_temperature(hydro_nt; output_path = nothing) !== nothing
        @test TethysDrift.plot_hydrodynamic_ts_diagram(hydro_nt; output_path = nothing) !== nothing

        trajs = (
            lons = [-63.5, -63.0],
            lats = [44.5, 45.0],
            depths = [-15.0, -25.0],
            stages = [:zoea_1, :settled_juvenile],
            alive = [true, true],
            settlement_status = [:active, :settled_successful],
            times = [0.0, 3600.0]
        )
        strata = [
            (name = "S1", lons = [-64.0, -63.0, -63.0, -64.0], lats = [44.0, 44.0, 45.0, 45.0])
        ]
        @test TethysDrift.plot_particle_trajectories(
            trajs;
            strata = strata,
            rotation_deg = 35.0,
            output_path = nothing
        ) !== nothing
    end

    @testset "Larval rules, tracking & storage round-trip" begin
        # Quadratic thermal mortality matching TethysBio
        @test isapprox(
            TethysDrift.larval_thermal_mortality_rate(11.0; base_mortality = 0.02),
            0.02;
            atol = 1e-12
        )
        @test isapprox(
            TethysDrift.larval_thermal_mortality_rate(5.0; base_mortality = 0.02),
            0.02 * (1.0 + 0.5 * (6.0 / 3.5)^2);
            atol = 1e-12
        )
        @test TethysDrift.larval_baseline_hazard(
            :instar1_settled;
            hazards = (0.020, 0.020, 0.015)
        ) == 0.0015
        @test TethysDrift.min_stage_days_for(
            :zoea1,
            TethysDrift.PELAGIC_INTERMOLT_DAYS
        ) == 36.0

        # Arrhenius-weighted Jan-Apr sole-cohort probabilistic hatching times
        @test isapprox(TethysDrift.arrhenius_hatch_rate(3.0), 1.0; atol = 1e-12)
        @test TethysDrift.arrhenius_hatch_rate(6.0) > 1.0
        @test 0.0 < TethysDrift.arrhenius_hatch_rate(0.0) < 1.0
        rt_warm = TethysDrift.sample_brood_hatch_times(
            fill(-62.0, 200), fill(44.5, 200), fill(-120.0, 200);
            temperature_fn = (lo, la, z, t) -> 6.0,
            origin_doy = 153.0,
            total_duration = 547.5 * 86400.0,
            rng = TethysDrift.Random.MersenneTwister(42)
        )
        rt_cold = TethysDrift.sample_brood_hatch_times(
            fill(-62.0, 200), fill(44.5, 200), fill(-120.0, 200);
            temperature_fn = (lo, la, z, t) -> 0.5,
            origin_doy = 153.0,
            total_duration = 547.5 * 86400.0,
            rng = TethysDrift.Random.MersenneTwister(42)
        )
        @test length(rt_warm) == 200
        @test (maximum(rt_warm) - minimum(rt_warm)) / 86400.0 <= 120.0
        @test (maximum(rt_cold) - minimum(rt_cold)) / 86400.0 <= 120.0
        @test (sum(rt_warm) / length(rt_warm)) < (sum(rt_cold) / length(rt_cold))  # Warmer benthic sites hatch earlier in Jan-Apr
        @test maximum(rt_cold) <= (547.5 - 150.0) * 86400.0 + 1.0

        # CFA polygons from docs/*.dat are loaded automatically in shelf order (North -> South -> 4X)
        cfa_polys = TethysDrift.load_cfa_polygons("nonexistent_inputs_dir")
        @test length(cfa_polys) == 3
        @test [p.code for p in cfa_polys] == [:cfanorth, :cfasouth, :cfa4x]
        @test length(cfa_polys[1].lons) == 15
        @test length(cfa_polys[2].lons) == 26
        @test length(cfa_polys[3].lons) == 28
        # Verify point-in-polygon connectivity uses docs/*.dat polygons
        conn_trajs = (
            lons = [-59.5 -59.5; -60.5 -60.5; -65.0 -65.0],
            lats = [46.5 46.5; 44.0 44.0; 42.5 42.5],
            depths = [-100.0 -100.0; -100.0 -100.0; -100.0 -100.0],
            stages = [:zoea1 :instar1_settled; :zoea1 :instar1_settled; :zoea1 :instar1_settled],
            alive = [true, true, true],
            survival_probability = [1.0 1.0; 1.0 1.0; 1.0 1.0],
            times = [0.0, 3600.0]
        )
        conn_res = TethysDrift.compute_empirical_connectivity(conn_trajs)
        @test conn_res.strata_names[1:3] == ["CFA North (CFA 20-22)", "CFA South (CFA 23-24)", "CFA 4X (Southwest NS)"]
        @test conn_res.counts_matrix[1, 1] == 1
        @test conn_res.counts_matrix[2, 2] == 1
        @test conn_res.counts_matrix[3, 3] == 1

        # End-to-end initialize_larval_particles + track_larval_cohort smoke test
        particles = TethysDrift.initialize_larval_particles(
            2;
            lon_range = (-63.2, -62.8),
            lat_range = (44.4, 44.6)
        )
        vel_fn = (lon, lat, z, t) -> (0.05, 0.0, 0.0)
        temp_fn = (lon, lat, z, t) -> 5.0
        bathy_fn = (lon, lat) -> -100.0
        res = TethysDrift.track_larval_cohort(
            particles;
            velocity_fn = vel_fn,
            temperature_fn = temp_fn,
            bathymetry_fn = bathy_fn,
            total_duration = 3600.0 * 6,
            dt = 1800.0,
            enable_molting = true,
            cv_mortality = 0.2
        )
        @test size(res.lons, 1) == 2
        @test size(res.stages, 1) == 2

        # Storage round-trip (open_storage, save_simulation_run!, list_simulation_runs, load)
        mktempdir() do tmp
            cfg = TethysDrift.load_configuration(
                joinpath(@__DIR__, "..", "configs", "default.toml")
            )
            opts = TethysDrift.configuration_to_options(cfg; config_name = "storage_test")
            db = TethysDrift.open_storage(joinpath(tmp, "db.zarr"))
            TethysDrift.save_simulation_run!(db, "run_1", opts; trajectories = res)
            @test db["trajectories"]["run_1"]["lons"].metadata.compressor isa Zarr.BloscCompressor
            @test db["trajectories"]["run_1"]["lons"].metadata.compressor.cname == "zstd"
            @test db["trajectories"]["run_1"]["lons"].metadata.compressor.clevel == 7
            @test db["trajectories"]["run_1"]["lons"].metadata.compressor.shuffle == 1
            @test eltype(db["trajectories"]["run_1"]["lons"]) == Float32
            runs_df = TethysDrift.list_simulation_runs(db)
            @test size(runs_df, 1) == 1
            loaded = TethysDrift.load_trajectories_namedtuple(db, "run_1")
            @test size(loaded.lons) == size(res.lons)
            TethysDrift.close_storage(db)

            html_path = joinpath(tmp, "tracks.html")
            sidecar_path = joinpath(tmp, "tracks.data.zz")
            TethysDrift.export_interactive_tracks_html(html_path; trajectories = res)
            @test isfile(html_path)
            @test isfile(sidecar_path)
            html_text = read(html_path, String)
            @test occursin("const SIDECAR_FILE = \"tracks.data.zz\";", html_text)
            @test occursin("id=\"nearest-n-input\"", html_text)
            @test occursin("function selectNearestReleases", html_text)
            @test occursin("id=\"ui-tooltip\"", html_text)
            @test occursin("data-tip=", html_text)
            @test occursin("window.location.protocol === \"file:\"", html_text)
            @test occursin("AbortController", html_text)
            @test occursin("rows.join(\"\\n\")", html_text)
            @test occursin("id=\"header-date\"", html_text)
            @test occursin("function formatNominalDate", html_text)
            @test !occursin("const RUNS = {", html_text)
            @test !occursin("Demographics & Zones", html_text)
            @test !occursin("const STRATA =", html_text)
            compressed = read(sidecar_path)
            @test length(compressed) >= 2 && compressed[1] == 0x78
            raw_json = String(Zarr.zuncompress(compressed, Zarr.ZlibCompressor(7), UInt8))
            parsed_runs = TethysDrift.JSON3.read(raw_json)
            first_run = first(values(parsed_runs))
            @test !haskey(first_run, :hydrodynamics)
            @test !haskey(first_run, :stats)
            @test !haskey(first_run, :connectivity_flows)
            @test length(first_run.particles) == size(res.lons, 1)
            @test length(first_run.times) == size(res.lons, 2)

            nc_path = joinpath(tmp, "out.nc")
            TethysDrift.export_larval_dispersal_netcdf(
                nc_path;
                trajectories = res,
                lon_bins = range(-64.0, -62.0; length = 4),
                lat_bins = range(44.0, 45.0; length = 4),
                config = cfg
            )
            @test isfile(nc_path)

            zarr_analysis_path = joinpath(tmp, "analysis.zarr")
            TethysDrift.export_larval_dispersal_jld2(
                zarr_analysis_path;
                trajectories = res,
                lon_bins = range(-64.0, -62.0; length = 4),
                lat_bins = range(44.0, 45.0; length = 4),
                config = cfg
            )
            @test isdir(zarr_analysis_path)

            hydro_snap = (
                lons = collect(-65.0:0.5:-62.0),
                lats = collect(43.0:0.5:45.0),
                depths = [-10.0, -50.0],
                u = fill(0.1, 7, 5, 2),
                v = fill(0.05, 7, 5, 2),
                speed = fill(0.11, 7, 5, 2),
                time_seconds = 0.0
            )
            hydro_snap.u[1, 1, 1] = NaN
            hydro_snap.v[1, 1, 1] = NaN
            hydro_snap.speed[1, 1, 1] = NaN
            adv_path = joinpath(tmp, "advection.png")
            TethysDrift.plot_hydrodynamic_advection(hydro_snap; output_path = adv_path)
            @test isfile(adv_path)
        end
    end

    @testset "Hydrodynamic time-origin alignment & snowcrab.toml" begin
        snow_cfg = TethysDrift.configuration_to_options(
            TethysDrift.load_config(joinpath(@__DIR__, "..", "configs", "snowcrab.toml"))
        )
        @test snow_cfg.track_duration ≈ 340.0 * 86400.0
        @test snow_cfg.hydro_start_days == 207.0

        # Domain-boundary exit rule marks exiting larvae as :off_domain and stops tracking them
        p_exit = TethysDrift.initialize_larval_particles(
            1; lon_range = (-63.0, -63.0), lat_range = (44.0, 44.0)
        )
        res_exit = TethysDrift.track_larval_cohort(
            p_exit;
            velocity_fn = (lo, la, z, t) -> (0.5, 0.0, 0.0),
            temperature_fn = (lo, la, z, t) -> 5.0,
            bathymetry_fn = (lo, la) -> -100.0,
            domain_interior_fn = (lo, la) -> lo <= -62.95,
            total_duration = 3600.0 * 12,
            dt = 1800.0,
            κ_h = 0.0,
            κ_v = 0.0,
            enable_molting = false
        )
        @test res_exit.settlement_status[1] == :off_domain
        @test res_exit.stages[1, end] == :off_domain
        @test res_exit.alive[1] == false

        mktempdir() do tmp
            zpath = joinpath(tmp, "offset_hydro.zarr")
            lons = [-64.0, -63.0]
            lats = [44.0, 45.0]
            depths = [-10.0, -50.0]
            # Simulated post-spinup clock starting at day 7 (604800 s), daily for 3 days
            times = [604800.0, 691200.0, 777600.0, 864000.0]
            root = Zarr.zgroup(zpath)
            _zput(g, name, arr) = (za = Zarr.zcreate(eltype(arr), g, name, size(arr)...); za[:] = arr)
            _zput(root, "lon", lons)
            _zput(root, "lat", lats)
            _zput(root, "depth", depths)
            _zput(root, "time", times)
            u = zeros(Float32, 2, 2, 2, 4)
            for m in 1:4
                u[:, :, :, m] .= Float32(m)
            end
            _zput(root, "u", u)
            _zput(root, "v", zeros(Float32, 2, 2, 2, 4))
            _zput(root, "w", zeros(Float32, 2, 2, 2, 4))
            _zput(root, "T", fill(5.0f0, 2, 2, 2, 4))

            store = TethysDrift.open_hydrodynamics(zpath)
            span = [0.0, 0.0]
            interp = TethysDrift.build_flow_from_store(
                store; start_time = times[1] + 86400.0, span = span
            )
            @test span[1] == 691200.0
            @test span[2] == 864000.0
            @test interp(-63.5, 44.5, -25.0, span[1]).u ≈ 2.0
            @test interp(-63.5, 44.5, -25.0, span[1] + 43200.0).u ≈ 2.5
        end
    end

    @testset "Wet-domain mask at spawn and advection" begin
        rng = TethysDrift.Random.MersenneTwister(123)
        wet_fn = (lon, lat, z) -> (lon >= -63.0 && lon <= -61.0 && lat >= 43.0 && lat <= 45.0)
        bathy_fn = (lon, lat) -> -120.0
        larvae = TethysDrift.initialize_larval_particles(
            20;
            lon_range = (-65.0, -60.0),
            lat_range = (42.0, 46.0),
            min_seabed_depth = 50.0,
            bathymetry = bathy_fn,
            wet_mask_fn = wet_fn,
            rng = rng
        )
        @test all(i -> wet_fn(larvae.lon[i], larvae.lat[i], larvae.depth[i]), 1:20)

        vel_fn = (lon, lat, z, t) -> wet_fn(lon, lat, z) ? (0.5, 0.5, 0.0) : (NaN, NaN, NaN)
        tr = TethysDrift.track_larval_cohort(
            larvae;
            velocity_fn = vel_fn,
            temperature_fn = (lon, lat, z, t) -> wet_fn(lon, lat, z) ? 5.0 : NaN,
            bathymetry_fn = bathy_fn,
            wet_mask_fn = wet_fn,
            total_duration = 3600.0 * 6,
            dt = 600.0,
            κ_h = 20.0,
            rng = rng
        )
        @test all(isfinite, tr.lons)
        @test all(isfinite, tr.lats)
        @test all(
            i -> wet_fn(tr.lons[i], tr.lats[i], tr.depths[i]),
            eachindex(tr.lons)
        )

        # Verify output_interval subsampling integrates at dt while recording at output_interval
        tr_sub = TethysDrift.track_larval_cohort(
            larvae;
            velocity_fn = vel_fn,
            temperature_fn = (lon, lat, z, t) -> wet_fn(lon, lat, z) ? 5.0 : NaN,
            bathymetry_fn = bathy_fn,
            wet_mask_fn = wet_fn,
            total_duration = 3600.0 * 6,
            dt = 600.0,
            output_interval = 3600.0 * 2,
            κ_h = 20.0,
            rng = TethysDrift.Random.MersenneTwister(123)
        )
        @test size(tr_sub.lons, 2) == 4
        @test tr_sub.times[1] == 0.0
        @test tr_sub.times[end] == 3600.0 * 6
        @test all(isfinite, tr_sub.lons)
    end
end
