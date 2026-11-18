# ==============================================================================
# Configuration, Domain Constants & Round-Trip Conversion Tests
# ==============================================================================

using Test
using TethysDrift
using Tethys

@testset "Configuration & Conversion" begin
    @testset "Domain config constants & discovery" begin
        @test TethysDrift.OCEAN_CONFIG == "ocean.toml"
        @test TethysDrift.BIOLOGY_CONFIG == "biology.toml"

        mktempdir() do tmp
            ocean_path = joinpath(tmp, TethysDrift.OCEAN_CONFIG)
            bio_path = joinpath(tmp, TethysDrift.BIOLOGY_CONFIG)
            write(ocean_path, "[domain]\nlon_min = -65.0\n")
            write(bio_path, "[biology]\nn_particles = 10\n")

            found = TethysDrift.find_domain_config_paths(tmp)
            @test abspath(found.ocean) == abspath(ocean_path)
            @test abspath(found.biology) == abspath(bio_path)
        end
    end

    @testset "Configuration <-> Options round-trip" begin
        cfg = TethysDrift.load_configuration(
            joinpath(@__DIR__, "..", "configs", "default.toml")
        )
        atmo = get!(cfg, "atmosphere", Dict{String, Any}())
        atmo["bulk_heat_flux"] = true
        atmo["source"] = "era5"
        stor = get!(cfg, "storage", Dict{String, Any}())
        delete!(stor, "ocean_model")
        hydro_sec = get!(cfg, "hydrodynamics", Dict{String, Any}())
        delete!(hydro_sec, "hydro_model_file")
        stor["storage_backend"] = "jld2"
        dvm = get!(cfg, "dvm", Dict{String, Any}())
        dvm["zoea1_day_depth"] = 45.0
        dvm["zoea1_night_depth"] = 8.0
        dvm["zoea2_depth_factor"] = 0.85
        dvm["megalopa_swim_factor"] = 1.0

        opts = TethysDrift.configuration_to_options(cfg; config_name = "audit_test")
        @test opts.bulk_heat_flux == true
        @test opts.atmospheric_source == :era5
        @test opts.storage_backend == "jld2"
        @test opts.zoea1_day_depth == 45.0
        @test opts.zoea1_night_depth == 8.0
        @test opts.zoea2_depth_factor == 0.85
        @test opts.megalopa_swim_factor == 1.0
        @test opts.hydro_model_file == joinpath(
            Tethys.lakehouse_dir("TethysFlow", "audit_test"),
            "hydrodynamics.zarr"
        )

        roundtrip = TethysDrift.options_to_configuration(opts)
        @test roundtrip["atmosphere"]["bulk_heat_flux"] == true
        @test roundtrip["storage"]["storage_backend"] == "jld2"
        @test roundtrip["dvm"]["zoea1_day_depth"] == 45.0
        @test roundtrip["dvm"]["zoea1_night_depth"] == 8.0
        @test roundtrip["dvm"]["zoea2_depth_factor"] == 0.85
        @test roundtrip["dvm"]["megalopa_swim_factor"] == 1.0
        @test !haskey(roundtrip["molting_and_settlement"], "cv_molt")
        @test !haskey(roundtrip["molting_and_settlement"], "cv_mortality")
    end

    @testset "to_hydrodynamic_config & to_larval_config" begin
        cfg = TethysDrift.load_configuration(
            joinpath(@__DIR__, "..", "configs", "default.toml")
        )
        opts = TethysDrift.configuration_to_options(cfg; config_name = "audit_test")
        hcfg = TethysDrift.to_hydrodynamic_config(opts)
        @test hcfg.nz_above == opts.nz_above
        @test hcfg.vertical_break_depth == opts.vertical_break_depth
        @test hcfg.cd_drag == 0.0025
        @test hcfg.storage_backend == opts.storage_backend
        @test hcfg.output_schedule_seconds == opts.output_schedule_seconds
        @test hcfg.checkpoint_schedule isa Float64

        lcfg = TethysDrift.to_larval_config(opts)
        @test lcfg.enable_zarr == opts.enable_zarr
        @test lcfg.zarr_path == opts.zarr_path
        @test lcfg.run_id == opts.run_id
        @test lcfg.seed == opts.seed
    end
end
