# Review Notes

## Main entry point

The 10k stochastic simulation entry point is the extensionless file `10k SAimulations`.
It defines:

```matlab
function Results = Overallfinding3(N_runs)
```

When called without arguments, it sets `N_runs = 10000`.

Because MATLAB expects a function file named after the first function, a valid-name copy was added as `Overallfinding3.m`. In the second pass, both `Overallfinding3.m` and the original extensionless file were updated to use the same shared deterministic metric helper.

## What the current 10k code does

For each Monte Carlo run, `Overallfinding3`:

1. Samples six kinetic rates from the existing log-uniform ranges:
   - `k1  = 10^( 3 + 4*rand)`
   - `k_1 = 10^(-1 + 4*rand)`
   - `k2  = 10^( 3 + 4*rand)`
   - `k_2 = 10^(-1 + 4*rand)`
   - `k3  = 10^(-1 + 4*rand)`
   - `k_3 = 10^( 1 + 4*rand)`
2. Sweeps 500 total-enzyme values from `10e-9` to `5.0e-6`.
3. Runs an `ode15s` solve at each enzyme value to compute final-time `P/(P+S)`.
4. Defines `E_star` as the enzyme value closest to the midpoint between the minimum and maximum ratio over that enzyme sweep.
5. Defines `E0 = 0.9 * E_star` and `E1 = 1.1 * E_star`, clipped to the enzyme sweep bounds.
6. Runs 200 Gillespie/SSA trajectories at `E0` and 200 at `E1`.
7. Computes KL divergence between the two final-time stochastic distributions using 40 histogram edges, PDF normalization, conversion back to bin probabilities, and epsilon smoothing with `eps_val = 1e-15`.
8. Computes energy dissipation, delta mu, and net flux at `E0` and `E1`.
9. Saves the full `Results` struct to `MonteCarloResults.mat`.

The diagnostic script `diagnose_Overallfinding3.m` keeps these definitions intact, except it defaults to `N_runs = 5` so timings can be collected before refactoring.

## Bottlenecks to measure

The dominant expected costs are:

- The deterministic midpoint search: 500 `ode15s` solves per parameter set.
- The stochastic sampling: 400 Gillespie trajectories per parameter set.
- Repeated steady-state ODE solves after the midpoint search:
  - `compute_entropy_local` solves ODE once for `E0` and once for `E1`.
  - `compute_dmu_local` solves ODE once for `E0` and once for `E1`.
  - `compute_Jnet_local` solves ODE once for `E0` and once for `E1`.

Before the second-pass refactor, each parameter set performed 506 ODE solves plus 400 stochastic trajectories.

The new diagnostic CSV records separate timings for:

- parameter sampling
- ODE sweep for `E_star`
- Gillespie simulations at `E0`
- Gillespie simulations at `E1`
- KL calculation
- energy calculation
- delta-mu calculation
- net-flux calculation
- combined energy/delta-mu/net-flux calculation
- saving/output
- total run time

## Scientific-definition concerns

These are notes only. No scientific definitions were changed.

- `E_star`, `E0`, and `E1` are currently defined by varying total enzyme concentration, not by varying substrate/input `S`.
- The paper discusses information rate and energy per burst cycle. This repo currently computes final-time KL divergence between stochastic output distributions and compares it with raw energy dissipation-like quantities.
- KL divergence depends on histogram binning and epsilon smoothing. The current 10k code uses 40 histogram edges, `eps_val = 1e-15`, and 200 stochastic samples per condition. Apparent KL saturation or sharp trends may be affected by binning and finite-sample artifacts.
- The stochastic simulator returns only final-time `P/(P+S)`. It does not estimate a time-resolved information rate.
- Energy dissipation, delta mu, and net flux are each computed from separate final-time ODE solves even when they use the same `Etot` and parameters.

## Files changed

- Added `Overallfinding3.m`: valid MATLAB filename copy of `10k SAimulations`.
- Updated `10k SAimulations`: original extensionless entry-point text now uses the shared metric helper for the post-KL deterministic metrics.
- Updated `Overallfinding3.m`: valid MATLAB entry point now uses the shared metric helper for the post-KL deterministic metrics.
- Added and updated `diagnose_Overallfinding3.m`: standalone small-run timing diagnostic that writes `Overallfinding3_diagnostics.csv` by default.
- Added `compute_steady_state_metrics_local.m`: shared deterministic helper that solves the ODE once for a parameter set and `Etot`, then computes species values, `P/(P+S)`, phi, delta mu, Jnet, and forward/backward flux diagnostics from that same final state.
- Added `check_metric_refactor_equivalence.m`: deterministic equivalence check comparing the old six-solve metric path with the new shared helper over small random parameter cases.
- Added `REVIEW_NOTES.md`: this review and refactor note.

## Second-pass redundancy refactor

The post-KL deterministic metric calculations in `Overallfinding3.m` now call `compute_steady_state_metrics_local` once for `E0` and once for `E1`.

The removed runtime redundancy is:

- old path: phi at `E0`, phi at `E1`, delta mu at `E0`, delta mu at `E1`, Jnet at `E0`, Jnet at `E1` each solved the same ODE separately;
- new path: one ODE solve at `E0` and one ODE solve at `E1`, with all formulas evaluated from those final states.

No changes were made to:

- parameter ranges;
- Gillespie/SSA behavior;
- KL binning;
- `N_stoch`;
- the definition of `E_star`, `E0`, or `E1`;
- the phi, delta-mu, Jnet, or ratio formulas.

The estimated ODE solve count per parameter set is now 502 instead of 506:

- 500 ODE solves for the unchanged `E_star` enzyme sweep;
- 2 ODE solves for shared deterministic metrics at `E0` and `E1`;
- 400 unchanged Gillespie trajectories.

This is only a small total reduction because the 500-point midpoint sweep still dominates deterministic solves, and the SSA work is unchanged.

## Equivalence checks

`check_metric_refactor_equivalence.m` checks:

- phi / energy dissipation;
- delta mu;
- Jnet / net flux.

It uses 3-5 small random parameter cases by default (`N_cases = 5`) and reports absolute and relative differences between:

- the old helper path, which solves separately for each metric;
- the new shared helper, which solves once and computes all metrics from one final state.

The diagnostic script also records per-run `max_metric_abs_diff` and `max_metric_rel_diff`, plus legacy and shared metric timings.

## Fast Smoke Diagnostic

`diagnose_Overallfinding3.m` now accepts optional name-value arguments or a config struct for diagnostic-only smoke testing:

```matlab
diagnose_Overallfinding3("N_runs",1,"n_E_vals",20,"N_stoch",5,"t_end_stoch",10,"useParallel",false)
```

The defaults still preserve the current diagnostic behavior:

- `N_runs = 5`;
- `n_E_vals = 500`;
- `N_stoch = 200`;
- `t_end_stoch = 1500`;
- `useParallel = false`;
- `use_shared_metrics = true`.

Any run with reduced `N_runs`, `n_E_vals`, `N_stoch`, or `t_end_stoch` is labeled in the output table as `fast_diagnostic_only` and includes a warning that the output is for smoke testing only and is not scientifically valid.

Fast smoke tests can produce identical final stochastic ratios because `N_stoch` and `t_end_stoch` are intentionally tiny. In that degenerate case only, `diagnose_Overallfinding3.m` widens the KL histogram edges just enough to avoid a histogram-edge error and records `KL_degenerate_bins_guard = true`.

I could not execute these checks in the current shell because neither `matlab` nor `octave` is installed on `PATH`, and no MATLAB app was discoverable under `/Applications`. Once MATLAB is available, run:

```matlab
check_metric_refactor_equivalence(5)
diagnose_Overallfinding3(5)
diagnose_Overallfinding3("N_runs",1,"n_E_vals",20,"N_stoch",5,"t_end_stoch",10,"useParallel",false)
```

Expected numerical differences should be zero or at roundoff level, since the ODE system and formulas are unchanged. Non-roundoff differences would need investigation before further optimization.

## Proposed next refactors

These should be considered after reviewing diagnostic timing output:

1. Review `check_metric_refactor_equivalence` and diagnostic CSV differences in MATLAB before making larger changes.
2. Decide whether it is scientifically acceptable to reuse final states from the 500-point midpoint sweep when `E0` or `E1` coincide with swept values. This was not done in the second pass.
3. Keep the SSA path intact, but add diagnostics for Gillespie step counts, early `max_steps` termination, and final-time molecule counts.
4. Add KL sensitivity diagnostics that run alongside the current default, such as alternate bin counts and bootstrap confidence intervals, without changing the default KL calculation.
5. Add explicit metadata to output files for parameter ranges, `N_stoch`, bin count, epsilon, `t_end`, `t_end_stoch`, and volume factor.
6. Only after the diagnostics are reviewed, consider parallel layout, caching, or vectorization changes with regression checks against the current implementation.
