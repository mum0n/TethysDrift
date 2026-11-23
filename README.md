# To see the simulation go here: [https://mum0n.github.io/TethysDrift/](https://mum0n.github.io/TethysDrift/)



# TethysDrift

Snow-crab larval-connectivity modelling. It tracks snow-crab larvae through a realized ocean
simulation and derives recruitment/connectivity metrics.

## What it does

`TethysDriftRun.jl` dispatches (`src/driver.jl`, library code) to:

- **`run_segment_tracking`** — release a larval cohort into an ocean store, advect/disperse/DVM/
  degree-day-molt the particles through it, and store trajectories.
- **`run_segment_metrics`** — recruitment/connectivity from the trajectories.
- **`render_all_figures`** — figures.

### What it reads

The advective field (`src/data/store_reader.jl`, `create_flow_interpolator`) opens a
hydrodynamic store and reads `u,v,w,T` per interval with the store's own cell-centre
coordinates. The spawn domain is taken from the store's own axes when present, CFA boundaries
otherwise.

That store is **input data produced elsewhere**: a `hydrodynamics.zarr` Tethys store written by
an ocean run — [`TethysFlow`](../TethysFlow/README.md) produces the ones this package's
configs default to. It is read-only here. The store's own coordinates are used rather than a
grid declared in this package, so the two cannot drift apart.

## Running

```
julia --project=. TethysDriftRun.jl --rerun \
  --hydro-file=<ocean store> \
  --config=default            # or --config=snowcrab
```

`--rerun` **or** `--force` run tracking; without a `--hydro-file` it needs one — this is the one
coupling, and it is a path, not a shared config. `--verbose`/`--viz` render figures from the
stored results.

## Configurations

Two, in `configs/` (the biology schema sections — `biology`, `dvm`,
`molting_and_settlement` — are the larval model's own):

| config | purpose |
|---|---|
| `default.toml` | Cape Breton demonstration — small cohort, generic placeholder biology, tracked through the keyless ocean `default.toml` store. |
| `snowcrab.toml` | Snow Crab production — 500-particle sole annual cohort over the 547.5-day (1.5-yr) `TethysFlow` span, with Arrhenius-weighted Jan–Apr probabilistic hatching (`50–300 m` shelf depth), `TethysBio`-aligned Poisson-CDF molting, DVM, and benthic settlement filtering. |

- **`[storage]`** is where the run reads and writes: `output_dir` (this package's results
  directory) and `ocean_model` (the realized ocean store to track through — the
  `--hydro-file` default). Both default into `lakehouse_dir("Tethys{Drift,Flow}", …)`;
  a command-line `--hydro-file` overrides. `ocean_model` names the Tethys store itself —
  the `hydrodynamics.zarr` directory the ocean half wrote — not the run directory holding
  it.
- A configuration that names no output directory of its own still writes into the
  lakehouse, under its own name. A run never lands in the package root: that is outside
  the lakehouse, where the other half never looks and no later run finds it.
- No `--config` ⇒ `configs/default.toml`.

### Shared Lakehouse Configuration

`TethysDrift` resolves all lakehouse directories via `Tethys.lakehouse_dir(...)`.
Precedence hierarchy:
- `ENV["TETHYS_LAKEHOUSE"]` or `ENV["GEODATA_LAKEHOUSE"]` (global override).
- `LocalPreferences.toml` (`[Tethys] lakehouse_dir`).
- `~/.tethys/storage.toml` or `~/.geodata/storage.toml` (`[storage] store_root`).
- Proximity fallback (`projects/lakehouse/`).

To configure the lakehouse persistently:
```julia
using Tethys
Tethys.set_lakehouse_dir!("/path/to/shared/lakehouse")
```

## Results

Trajectories, metrics, shelf-aligned (35° rotated) PNG diagnostic maps, and the canvas-backed interactive Leaflet dashboard (`interactive_larval_tracks.html` + zlib sidecar `.data.zz`) are written to `[storage] output_dir` — separate from the read-only input ocean store.

## Documentation

- [`docs/TethysDrift.md`](docs/TethysDrift.md) — full reference manual (biophysical equations, configuration, analysis, storage, visualization).
- [`docs/TethysDriftNotes.md`](docs/TethysDriftNotes.md) — concise biophysical equations, hatching/molting rules, and assumptions.
- [`docs/github_pages.md`](docs/github_pages.md) — local viewing and GitHub Pages deployment for the interactive Leaflet dashboard.
- [`docs/snow_crab_larval_connectivity_paper.md`](docs/snow_crab_larval_connectivity_paper.md) — snow-crab larval connectivity manuscript.
