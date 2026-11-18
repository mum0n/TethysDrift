# TethysDrift: Biophysical Larval Transport and Demographic Connectivity

Technical reference for the Lagrangian larval drift and demographic connectivity workflow for snow crab (*Chionoecetes opilio*) across the Scotian Shelf and Gulf of St. Lawrence.

## Pipeline Architecture and Scientific Workflow

```text
[TethysFlow Hydrodynamic Store] (u, v, w, T, S from Zarr)
               │
               ▼
   [Stage 1: Cohort Initialization] ── Benthic hatching along Crab Fishing Areas (CFAs)
               │
               ▼
   [Stage 2: Active Vertical Ascent] ── Upward swimming (10 mm/s) to surface
               │
               ▼
   [Stage 3: Pelagic Advection & SDE] ── Ocean currents + Stokes drift + Horizontal/Vertical diffusion
               │
               ├── Context-dependent Diel Vertical Migration (DVM)
               ├── Temperature-dependent degree-day accumulation (DD)
               └── Multi-stage ecdysis (Zoea I -> Zoea II -> Megalopa)
               │
               ▼
   [Stage 4: Benthic Nursery Settlement] ── Substrate depth, temperature ceiling (< 6°C), tidal filtering
               │
               ▼
   [Stage 5: Demographic Connectivity] ── Transition matrices, self-recruitment, source-sink dynamics
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
- **Thermal Base $T_0 = -1.5^\circ\text{C}$**: Physiological development halts at $-1.5^\circ\text{C}$; below this threshold, larvae undergo thermal quiescence.
- **Benthic Settlement Thermal Window**: Settlement is permitted only between $-350\text{ m} \le z \le -50\text{ m}$ where bottom temperatures satisfy $T_{\text{bed}} \le 6.0^\circ\text{C}$. Energetic tidal fluctuations are filtered using a 12.42-hour moving exponential average.
- **Multi-Scale Spatio-Temporal Coupling**:
  - **Ocean Hydrodynamics (`TethysFlow`)**: Solves 3D Navier-Stokes equations with barotropic subcycling at $\Delta t = 30 - 60\text{ s}$ on a structured $\sim 2.5\text{ km}$ spherical grid, resolving the baroclinic Rossby radius ($R_d \approx 4 - 8\text{ km}$). Output is archived at a 1-hour cadence ($\Delta t_{\text{output}} = 3600\text{ s}$) to capture semidiurnal tidal excursions ($M_2 \approx 12.42\text{ h}$) and inertial ellipses ($T_f \approx 17.1\text{ h}$).
  - **Offline Lagrangian Advection (`TethysDrift`)**: Particles are stepped via RK4 at $\Delta t_{\text{part}} = 300\text{ s}$ ($5\text{ min}$) by spatiotemporal interpolation of the hourly flow tensor, supplemented by subgrid turbulent dispersion ($K_h \approx 10\text{ m}^2/\text{s}$). Bathymetric steering ($\sim 100\text{--}500\text{ m}$) is enforced locally at the seabed boundary.

### Management and Ecological Utility

- **Marine Protected Area (MPA) Network Design**: Identifies persistent larval source zones (e.g., St. Anns Bank, Gully MPA) and recipient nurseries.
- **Stock Assessment and Spatially Explicit Recruitment**: Provides empirical recruitment indices to calibrate population dynamics models (`TethysBio`).
- **Climate Vulnerability Projections**: Quantifies how ocean warming contracts pelagic larval duration (PLD) and alters downstream transport pathways.
