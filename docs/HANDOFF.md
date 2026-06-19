# Handoff

Short deliverable summary for the next person (researcher / supervisor). For the
full reasoning see [MODEL_AUDIT.md](MODEL_AUDIT.md) and
[RESULT_INTERPRETATION.md](RESULT_INTERPRETATION.md).

## What was fixed

- **KL ~ 34 saturation artifact removed.** The old plateau was a numerical
  artifact of histogram KL with `epsilon = 1e-15` (ceiling `-log(epsilon)`),
  made worse by asymmetric high-driving rate ranges and RT-scaled units. It is
  gone in the dimensionless runs.
- **Staged, reproducible pipeline** (screen -> select -> SSA/KL -> plot) replaces
  the ~48 h monolithic 10k brute force. Runtime for a real dimensionless result
  is now ~1-2 minutes.
- **`paper_inspired_dimensionless` mode**: ATP=1, RT=1 (phi/delta_mu in kBT-like
  units), symmetric log-uniform rate ranges 10^-2..10^2, E0=0.95*E*, E1=1.05*E*.
- **Bounded, interpretable metrics**: raw KL kept for continuity, plus
  `symmetric_KL`, `log1p_symmetric_KL`, and **JS divergence (the headline,
  bounded 0..ln2)**.
- **Validity filtering and diagnostics**: ODE failures, complex/negative phi,
  `max_steps_hit`, sparse/degenerate histograms, NaN/Inf, plus boundary-fraction
  stats and Spearman correlations reported in `summary.md`.
- **Log-spaced enzyme sweep option** (`enzyme_spacing='log'`) that places E* on
  genuine response transitions (`E_star_at_sweep_boundary` -> 0).
- **Legacy mode preserved exactly** for comparison.

## What outputs were generated

| Location | Contents |
|----------|----------|
| `results/dimensionless_first_result/` | linear-sweep CSVs + `summary.md` |
| `figures/dimensionless_first_result/` | linear-sweep figures (12 PNGs) |
| `results/dimensionless_logspacing_first_result/` | log-sweep CSVs + `summary.md` (**researcher-facing**) |
| `figures/dimensionless_logspacing_first_result/` | log-sweep figures (12 PNGs) |
| `results/dimensionless_smoke/`, `figures/dimensionless_smoke/` | fast smoke-test outputs |
| `results/first_result/`, `results/smoke/` (+ figures) | legacy-mode outputs (comparison) |

Key figures (in each dimensionless figures dir):

- `fig_JS_vs_phi.png` — **main** distinguishability-vs-energy plot.
- `fig_log1p_symmetric_KL_vs_phi.png`, `fig_symmetric_KL_vs_average_phi.png`,
  `fig_KL_vs_phi.png` — KL variants for comparison.
- `fig_example_response_curves.png` — deterministic P/(P+S) vs total enzyme with
  E0/E*/E1 marked (semilogx for the log sweep).
- `fig_example_histograms_by_energy.png` — low/medium/high E0 vs E1 output dists.
- `fig_validity_summary.png` — valid vs excluded SSA rows.

## What remains scientifically unresolved

The numerical artifacts are fixed, but **the enzyme model does not reproduce a
strong monotonic paper-like energy/discrimination trend** (Spearman |rho| <= ~0.4,
mildly negative), and log spacing did not change that. The remaining mismatch is
most likely **conceptual**. Four decisions are needed from the researcher /
supervisor (also in MODEL_AUDIT.md section 8):

1. Should the input signal be total enzyme concentration, substrate
   concentration, substrate influx alpha, or something else?
2. Is the enzyme mass-action model intended as an analogy, or should the actual
   PNAS transcriptional Markov model be implemented?
3. Is the desired output metric KL/JS between stochastic product distributions,
   or the paper's information-rate/decision metric?
4. Should E0/E1 be +/-5% around the midpoint (delta_c = 0.1*c*), or something else?

Current code assumes: total enzyme, analogy, KL/JS of SSA outputs, +/-5%.

A focused diagnostic (`results/diagnostics/DIAGNOSTIC_REPORT.md`, plus
`figures/diagnostics/`) localized *why* the trend is flat:

- **JS is decoupled from the deterministic response geometry** — high-JS rows can
  be flat, the lowest-JS rows can be sharp sigmoids — so it is **not** an
  E*-placement / sweep-resolution problem.
- **Raw P and S are single-digit counts**, and a `Vol_factor = 1/10/100` sweep
  collapses the 0/1 boundary pile-up (median boundary fraction 0.74 -> 0.14 ->
  0.00). So the boundary-heavy histograms are a **low-copy-number artifact**, not
  the metric `P/(P+S)` itself.
- **Bottom line:** the flat trend is *partly* a Vol_factor/discreteness artifact
  (fixable) and *primarily* an experiment-design issue (the input is ±5% total
  enzyme with substrate influx fixed, so energy is not an independent control).

## Reproduce the dimensionless result

From the repository root (MATLAB on PATH):

```bash
# fast smoke test (minutes)
matlab -batch "run('scripts/run_dimensionless_smoke_pipeline.m')"

# linear enzyme sweep (~1-2 min)
matlab -batch "run('scripts/run_dimensionless_first_result_pipeline.m')"

# log-spaced enzyme sweep — researcher-facing (~1-2 min)
matlab -batch "run('scripts/run_dimensionless_logspacing_first_result_pipeline.m')"

# focused diagnostics (no large runs; need the log-spaced result above)
matlab -batch "run('scripts/run_dimensionless_volfactor_sensitivity.m')"
matlab -batch "addpath('scripts'); setup_repo_paths(); diagnose_dimensionless('ssaFile','results/dimensionless_logspacing_first_result/results_ssa_kl.csv','enzyme_spacing','log')"
```

Each script calls `setup_repo_paths()` itself, so no manual path setup is needed.
Use only rows with `valid_for_kl = true` for any KL/JS-vs-energy interpretation.
Do not raise `t_end_stoch` or `N_stoch` to legacy brute-force values without a
convergence check.
