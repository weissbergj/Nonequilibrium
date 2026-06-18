# Result interpretation guide

## Scientific caveat (read first)

**E0/E1 are currently defined by total enzyme concentration.** The code sweeps total enzyme E, finds E* at the midpoint of the deterministic output ratio P/(P+S), then sets E0 = 0.9·E* and E1 = 1.1·E*.

If the intended input signal is substrate S or influx alpha, the scientific definition should be changed before final interpretation.

## Output files

### `results_*.csv` (deterministic screen)

One row per sampled parameter set.

| Column | Meaning |
|--------|---------|
| `k1` … `k_minus3` | Sampled rate constants |
| `E_star`, `E0`, `E1` | Enzyme input states (total enzyme, M) |
| `phi_E0`, `phi_E1` | Energy dissipation-like quantity at E0/E1 |
| `delta_mu_E0`, `delta_mu_E1` | Cycle chemical driving (J/mol) |
| `Jnet_E0`, `Jnet_E1` | Net flux at E0/E1 |
| `average_phi`, `average_delta_mu` | Mean of E0 and E1 values |
| `valid_for_screen` | Row passed deterministic quality filters |
| `ode_tolerance_warning` | ode15s tolerance warning during sweep/metrics |
| `complex_phi` | phi has imaginary part (reverse flux dominated step) |
| `negative_phi` | phi < 0 |
| `E_star_at_sweep_boundary` | E* at edge of sweep — E0/E1 may be unreliable |

**Use `valid_for_screen=true` rows for selection.**

### `selected_parameter_sets.csv`

Subset chosen for SSA, with `selection_bin` = `low_energy`, `medium_energy`, or `high_energy`.

### `results_ssa_kl.csv`

Deterministic columns plus stochastic/KL results.

| Column | Meaning |
|--------|---------|
| `KL_E0_E1`, `KL_E1_E0` | Histogram KL divergences |
| `symmetric_KL` | 0.5·(KL_E0_E1 + KL_E1_E0) — primary y-axis for trend plots |
| `JS_divergence` | Jensen–Shannon (bounded 0–ln2); supplementary sanity check |
| `output_mean_E0/E1`, `output_var_E0/E1` | Mean/var of final P/(P+S) across SSA trajectories |
| `mean_steps_E0/E1` | Average Gillespie steps per trajectory |
| `max_steps_hits_E0/E1` | Count of trajectories hitting step cap |
| `valid_for_kl` | Row acceptable for KL trend plots |
| `ssa_warning_flags` | Semicolon-separated issues |

**`valid_for_kl=false` when any of:**

- max_steps hit on any trajectory
- failed trajectory
- degenerate or sparse histogram
- N_stoch < 5·nBins (insufficient samples per bin)
- NaN/Inf/complex KL

### `ssa_example_E0_E1_distributions.csv`

Final P/(P+S) values for all SSA trajectories at E0 and E1 for the **first selected row**. Used for `fig_example_E0_E1_histograms.png`.

## Figures

| Figure | Description |
|--------|-------------|
| `fig_KL_vs_phi.png` | KL(E0‖E1) vs average phi |
| `fig_KL_vs_delta_mu.png` | KL vs average delta mu |
| `fig_symmetric_KL_vs_average_phi.png` | Primary energy–distinguishability plot |
| `fig_energy_distribution.png` | Histogram of average phi |
| `fig_delta_mu_distribution.png` | Histogram of average delta mu |
| `fig_example_E0_E1_histograms.png` | Overlaid E0/E1 output distributions |
| `fig_validity_summary.png` | Counts of valid vs invalid SSA rows |

Blue points = `valid_for_kl=true`. Gray = invalid (diagnostic only).

## How to read KL-vs-energy plots

1. Check `summary.md` for how many rows are `valid_for_kl`.
2. If most KL values cluster at a high plateau (~30+), suspect histogram artifact — reduce `nBins` or increase `N_stoch`.
3. Compare `symmetric_KL` with `JS_divergence`; large disagreement suggests unstable binning.
4. Check `max_steps_hits` — if high, SSA did not reach `t_end_stoch`; KL is not a converged steady-state measurement.
5. Compare `output_mean_E0` and `output_mean_E1` — if nearly equal but KL is huge, distrust KL.

## Recommended settings (current)

| Setting | Smoke | First result | Legacy (avoid) |
|---------|-------|--------------|----------------|
| N_runs | 50 | 500 | 10000 |
| n_E_vals | 50 | 100 | 500 |
| N_select | 12 | 50 (+ oversample) | all |
| N_stoch | 50 | 50 | 200 |
| t_end_stoch | 25 | 25 | 1500 |
| nBins | 10 | 10 | 40 |

Increase `N_stoch` and `t_end_stoch` only after confirming SSA output distributions stabilize.

## What to report to collaborators

1. Number of valid deterministic and SSA rows
2. KL-vs-phi plot using valid rows only
3. Example E0/E1 histogram
4. Explicit note that input is total enzyme, not substrate
5. SSA settings used (N_stoch, t_end_stoch, nBins)
6. Count of excluded rows and why (from `summary.md` and `ssa_warning_flags`)
