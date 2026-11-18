using Test
using TethysDrift
const TethysDrift = TethysDrift

@testset "Driver Segments & Surface API" begin
    # These lived in the umbrella's script, where `Pkg.test` could not reach them. Being
    # defined in `src/` is what makes them testable at all.
    for name in ("run_segment_tracking", "run_segment_metrics", "run_segment_visualize",
                 "render_all_figures", "load_run_trajectories", "run_cli_list_runs",
                 "run_cli_compare_runs", "run_cli_model_average",
                 "create_flow_interpolator", "open_hydrodynamics")
        @test isdefined(TethysDrift, Symbol(name))
    end

    @testset "CLI parse_args" begin
        cli_mod = Module(:TethysDriftCLI)
        Base.include(cli_mod, joinpath(@__DIR__, "..", "TethysDriftRun.jl"))
        parsed = Base.invokelatest(
            cli_mod.parse_args,
            [
                "--config=.\\configs\\snowcrab.toml",
                "--hydro_path=..\\lakehouse\\TethysFlow\\scotian_shelf\\",
                "--all"
            ]
        )
        @test parsed["verb"] == "all"
        @test parsed["config"] == ".\\configs\\snowcrab.toml"
        @test parsed["hydro_file"] == "..\\lakehouse\\TethysFlow\\scotian_shelf\\"

        for flag in ("--hydro-file=a.zarr", "--hydro_file=a.zarr",
                     "--hydro-path=a.zarr", "--hydro_path=a.zarr")
            p2 = Base.invokelatest(cli_mod.parse_args, [flag, "--rerun"])
            @test p2["hydro_file"] == "a.zarr"
            @test p2["verb"] == "rerun"
        end

        # Without verb or hydro flag, parse_args leaves verb=nothing (main defaults to --all)
        # and hydro_file="" (main resolves from config `[storage] ocean_model`).
        p_default = Base.invokelatest(cli_mod.parse_args, ["--config=.\\configs\\snowcrab.toml"])
        @test p_default["verb"] === nothing
        @test p_default["hydro_file"] == ""

        cfg = TethysDrift.load_configuration(
            joinpath(@__DIR__, "..", "configs", "snowcrab.toml")
        )
        opts = TethysDrift.configuration_to_options(cfg; config_name = "snowcrab")
        resolved_hydro, _ = TethysDrift.resolve_hydro_model_path(opts, "hydrodynamics.zarr")
        @test endswith(
            replace(resolved_hydro, "\\" => "/"),
            "lakehouse/TethysFlow/scotian_shelf/hydrodynamics.zarr"
        )
    end
end
