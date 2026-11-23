# TethysDrift: Biophysical Larval Transport and Demographic Connectivity

Technical reference for the Lagrangian larval drift and demographic connectivity workflow for snow crab (*Chionoecetes opilio*) across the Scotian Shelf and Gulf of St. Lawrence.

## Pipeline Architecture and Scientific Workflow

```text
[TethysFlow Hydrodynamic Store] (u, v, w, T, S from Zarr)
               │
               ▼
   [Stage 1: Cohort Initialization] ── Shelf seabed seeding (50-300 m) & Arrhenius Jan-Apr hatching
               │
               ▼
   [Stage 2: Active Vertical Ascent] ── Upward swimming (10 mm/s) to surface
               │
               ▼
   [Stage 3: Pelagic Advection & SDE] ── Ocean currents + Horizontal/Vertical diffusion
               │
               ├── Stage-dependent Diel Vertical Migration (DVM)
               ├── Per-stage degree-day accumulation (270 DD, base -1.5 °C)
               └── Poisson-CDF ecdysis hazard (Zoea I -> Zoea II -> Megalopa, aligned with TethysBio)
               │
               ▼
   [Stage 4: Benthic Nursery Settlement] ── Depth (10-500 m), temperature ceiling (<= 7.0 °C), M2 tidal filter
               │
               ▼
   [Stage 5: Demographic Connectivity] ── Exact docs/cfa*.dat polygon transitions, self-retention, export
```

## Biophysical Governing Equations and Core Assumptions

The trajectory of each individual larva $p$ is integrated using the stochastic Langevin equation via the Euler-Maruyama scheme:

$$
\boldsymbol{x}_h^{n+1} = \boldsymbol{x}_h^n + \left[ f_{\text{bbl}}(z^n) \boldsymbol{u}_h(\boldsymbol{x}^n, t^n) \right] \Delta t + \sqrt{2 \kappa_h \Delta t} \, \boldsymbol{\xi}_h^n
$$

The horizontal velocity $\boldsymbol{u}_h$ is read from the ocean half's store, which already carries the tidal cycle driven by TPXO9 boundary forcing, so no separate tidal term is added in the Lagrangian tracker.

$$
z^{n+1} = z^n + \left[ w(\boldsymbol{x}^n, t^n) + w_{\text{swim}}(z^n, t^n) + w_{\text{sink}}(\text{stage}) + \frac{\partial \kappa_v}{\partial z} \right] \Delta t + \sqrt{2 \kappa_v \left(z^n + \frac{1}{2}\frac{\partial \kappa_v}{\partial z}\Delta t\right) \Delta t} \, \xi_z^n
$$

### Core Formulations

- **Bottom Boundary Layer (BBL) Logarithmic Shear**:
  $$
  f_{\text{bbl}}(z) = \text{clamp}\left( \frac{\ln(\max(z_0, z - z_{\text{bed}}) / z_0)}{\ln(h_{\text{bbl}} / z_0)}, 0.0, 1.0 \right)
  $$
  where $h_{\text{bbl}} = 10.0\text{ m}$ and roughness length $z_0 = 10^{-3}\text{ m}$.
- **Deterministic Visser Diffusive Drift**: The gradient term $\frac{\partial \kappa_v}{\partial z}$ prevents artificial accumulation of particles in low-diffusivity pycnoclines (Visser, 1997).
- **Passive Gravitational Settling**: Stage-dependent sinking rates are Zoea I ($-0.5\text{ mm/s}$), Zoea II ($-1.0\text{ mm/s}$), and Megalopa ($-2.5\text{ mm/s}$).

## Assumptions and Scientific Utility

### Explicit Assumptions

- **Passive Horizontal Advection with Turbulent Dispersion**: Horizontal swimming motility ($< 5\text{ mm/s}$) is small relative to regional currents ($0.10 - 0.50\text{ m/s}$), justifying passive horizontal advection.
- **Deterministic Visser Correction**: Subgrid vertical mixing requires $\frac{\partial \kappa_v}{\partial z}$ to satisfy mass conservation and avoid unphysical accumulation in pycnoclines.
- **Arrhenius-Weighted Sole-Cohort Hatching (Jan 1 – Apr 30)**: Eggs are seeded on the continental shelf (`50–300 m` depth) and held stationary in `:egg` stage until their individual hatch time, drawn from a truncated cumulative Normal hazard on DOY `1–120` whose modal day shifts earlier in warmer bottom water via $r(\bar{T}_{\text{bed}}) = \exp(-(E_a/k_B)(1/(\bar{T}_{\text{bed}}+273.15) - 1/(T_{\text{ref}}+273.15)))$ ($E_a = 0.55\text{ eV}$, $T_{\text{ref}} = 3.0\,^\circ\text{C}$, $Q_{10} \approx 2.34$), aligned with `TethysBio`.
- **Arrhenius Physiological Degree-Days & `TethysBio`-Aligned Molting**: Per-stage physiological degree-days accumulate via $\Delta DD = (T_{\text{ref}} - T_0) \cdot r(T) \, \Delta t$ anchored at $T_{\text{ref}} = 3.0^\circ\text{C}$ and $T_0 = -1.5^\circ\text{C}$ (reset at each molt, without benthic thermal inhibition) and drive stage transitions via the continuous Poisson-CDF daily hazard $p_{\text{trans}} = 1 - \exp(-h_0 \cdot \text{PoissonCDF}(270; \lfloor DD \rfloor) \cdot \Delta t)$ with $h_0 = 0.35\text{ day}^{-1}$.
- **Mortality**: Pelagic baseline natural mortality is configurable (disabled as `[0, 0, 0]` in `snowcrab.toml`); mortality occurs via settlement in unsuitable habitat or the pelagic over-age developmental penalty (`0.05` per 10 d past 84 d in a single stage).
- **Benthic Settlement Window**: In `snowcrab.toml`, settlement is permitted between $-500\text{ m} \le z_{\text{bed}} \le -10\text{ m}$ where M2-filtered bottom temperature satisfies $\bar{T}_{\text{bed}} \le 7.0^\circ\text{C}$ (12.42-hour exponential filter).
- **Administrative CFA Polygons**: Transition probability matrices $P_{ij}$ classify release and settlement coordinates by ray-casting against exact polygon rings in `docs/cfanorth.dat`, `docs/cfasouth.dat`, and `docs/cfa4x.dat`, with a catch-all `Outside / Offshore` stratum.
- **Multi-Scale Spatio-Temporal Coupling**:
  - **Ocean Hydrodynamics (`TethysFlow`)**: Solves 3D Navier-Stokes equations with barotropic subcycling at $\Delta t = 30 - 60\text{ s}$ on a structured $\sim 2.5\text{ km}$ spherical grid, resolving the baroclinic Rossby radius ($R_d \approx 4 - 8\text{ km}$). Output is archived at a 1-hour cadence ($\Delta t_{\text{output}} = 3600\text{ s}$) to capture semidiurnal tidal excursions ($M_2 \approx 12.42\text{ h}$) and inertial ellipses ($T_f \approx 17.1\text{ h}$).
  - **Offline Lagrangian Advection (`TethysDrift`)**: Particles are stepped via RK4 at $\Delta t_{\text{part}} = 300\text{ s}$ ($5\text{ min}$) by spatiotemporal interpolation of the hourly flow tensor, supplemented by subgrid turbulent dispersion ($K_h \approx 10\text{ m}^2/\text{s}$). Bathymetric steering ($\sim 100\text{--}500\text{ m}$) is enforced locally at the seabed boundary.

### Management and Ecological Utility

- **Marine Protected Area (MPA) Network Design**: Identifies persistent larval source zones (e.g., St. Anns Bank, Gully MPA) and recipient nurseries.
- **Stock Assessment and Spatially Explicit Recruitment**: Provides empirical recruitment indices to calibrate population dynamics models (`TethysBio`).
- **Climate Vulnerability Projections**: Quantifies how ocean warming contracts pelagic larval duration (PLD) and alters downstream transport pathways.
