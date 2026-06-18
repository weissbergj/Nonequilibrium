# Nonequilibrium enzyme pipeline

Staged, reproducible MATLAB pipeline for studying KL divergence between stochastic output distributions at two nearby enzyme input states (E0, E1), versus deterministic energy dissipation (phi) and chemical driving (delta mu).

**Important:** E0/E1 are currently defined by **total enzyme concentration**, not substrate S. See [docs/RESULT_INTERPRETATION.md](docs/RESULT_INTERPRETATION.md).

## Quick start

From the repository root (MATLAB R2019b+ recommended):

```matlab
run('scripts/run_smoke_pipeline.m')
```

Or non-interactive:

```bash
matlab -batch "run('scripts/run_smoke_pipeline.m')"
```

Smoke run completes in about **5–10 minutes** and writes:

- `results/smoke/*.csv`
- `figures/smoke/*.png`
- `results/smoke/summary.md`

## First real run (recommended)

```matlab
run('scripts/run_first_result_pipeline.m')
```

Bounded settings: 500 deterministic screens, 50 stratified SSA rows, `N_stoch=50`, `t_end_stoch=25`, `nBins=10`. Expect **~30–90 minutes** depending on hardware.

Outputs: `results/first_result/`, `figures/first_result/`, `results/first_result/summary.md`.

## Repository layout

| Path | Purpose |
|------|---------|
| `legacy/` | Original extensionless files and brute-force `Overallfinding3` |
| `src/core/` | ODE model, SSA, KL, metrics, parameter sampling |
| `src/pipeline/` | Screen → select → SSA/KL → plot |
| `src/diagnostics/` | Equivalence checks |
| `scripts/` | Entry-point runners |
| `results/` | Generated CSV outputs |
| `figures/` | Generated PNG plots |
| `docs/` | Pipeline notes and interpretation guide |

## What to send the researcher

After `run_first_result_pipeline.m`:

1. `results/first_result/results_ssa_kl.csv` — KL, phi, validity flags
2. `results/first_result/selected_parameter_sets.csv` — chosen parameter sets with E*, E0, E1
3. `figures/first_result/fig_symmetric_KL_vs_average_phi.png` (and related plots)
4. `figures/first_result/fig_example_E0_E1_histograms.png`
5. `results/first_result/summary.md`
6. `docs/RESULT_INTERPRETATION.md`

Use only rows with **`valid_for_kl = true`** for KL-vs-energy trends.

## Full-scale template

`scripts/run_full_pipeline_template.m` is **disabled by default**. Uncomment and edit after SSA feasibility is confirmed. Do not use legacy defaults (`t_end_stoch=1500`, `nBins=40`) without validation.

## Documentation

- [docs/PIPELINE_NOTES.md](docs/PIPELINE_NOTES.md) — staged workflow vs old 48h brute-force
- [docs/RESULT_INTERPRETATION.md](docs/RESULT_INTERPRETATION.md) — CSV columns, plots, caveats

## Legacy brute-force script

The original monolithic 10k run lives under `legacy/10k SAimulations` and `legacy/Overallfinding3.m`. It is preserved for reference but is **not** the recommended workflow.
