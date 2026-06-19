# Nonequilibrium enzyme pipeline

Staged, reproducible MATLAB pipeline for studying KL/JS divergence between stochastic output distributions at two nearby enzyme input states (E0, E1), versus deterministic energy dissipation (phi) and chemical driving (delta mu).

**Important:** E0/E1 are currently defined by **total enzyme concentration**, not substrate S. See [docs/RESULT_INTERPRETATION.md](docs/RESULT_INTERPRETATION.md).

## Status / headline conclusion

The original "KL saturates at ~34" result was a **numerical artifact** (histogram KL + epsilon smoothing). The refactor + `paper_inspired_dimensionless` mode + JS/log-KL metrics + symmetric rate ranges + log-spaced enzyme sweep **remove that artifact**. However, the corrected enzyme model **does not show a strong monotonic paper-like energy/discrimination trend**, and log-spacing the sweep did not fix that. The remaining mismatch is most likely **conceptual** (model choice, input definition, or metric definition), not numerical.

See [docs/HANDOFF.md](docs/HANDOFF.md) for the deliverable summary and [docs/MODEL_AUDIT.md](docs/MODEL_AUDIT.md) for the full audit and the open questions for the researcher/supervisor.

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

## Paper-inspired dimensionless mode (recommended for interpretation)

Same enzyme ODE/SSA model, but dimensionless and paper-inspired: `ATP=1`, `RT=1`
(phi/delta_mu in kBT-like units), symmetric log-uniform rate ranges 10^-2..10^2,
`E0=0.95*E*`, `E1=1.05*E*`, bounded SSA. The **headline metric is JS divergence**
(bounded), not raw KL.

```bash
# linear enzyme sweep (~1 min)
matlab -batch "run('scripts/run_dimensionless_first_result_pipeline.m')"

# log-spaced enzyme sweep — cleaner E* definition, researcher-facing (~1 min)
matlab -batch "run('scripts/run_dimensionless_logspacing_first_result_pipeline.m')"
```

Outputs land under `results/dimensionless_first_result/`,
`results/dimensionless_logspacing_first_result/`, and the matching `figures/`
dirs. The main figure is `fig_JS_vs_phi.png`; `fig_example_response_curves.png`
shows the deterministic response and the E0/E\*/E1 operating points.

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
