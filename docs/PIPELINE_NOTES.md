# Pipeline notes

## Why not the old 10k brute-force run?

The legacy script `legacy/10k SAimulations` / `legacy/Overallfinding3.m`:

1. Sampled 10,000 random kinetic parameter sets.
2. For each set: 500 ODE solves to find E*, then 400 Gillespie trajectories, then KL.
3. Took ~48 hours and produced unclear KL-vs-energy plots.

The expensive step is Gillespie/SSA. Running SSA on every screened parameter set before knowing whether deterministic metrics are usable wastes most of the runtime.

## Staged workflow

| Stage | Function | Cost |
|-------|----------|------|
| 1 | `run_deterministic_screen` | ODE only |
| 2 | `select_parameter_sets` | bookkeeping |
| 3 | `run_ssa_kl_on_selected` | SSA + KL (expensive) |
| 4 | `plot_pipeline_results` | plotting |

### Stage 1 — deterministic screen

Computes per parameter set:

- Sampled rates k1…k_3 (legacy log-uniform ranges, unchanged)
- E*, E0 = 0.9·E*, E1 = 1.1·E* (total enzyme sweep, unchanged)
- phi, delta mu, Jnet at E0 and E1
- Validity flags: `valid_for_screen`, `ode_tolerance_warning`, `complex_phi`, etc.

Does **not** compute KL (KL requires stochastic distributions).

### Stage 2 — selection

Excludes invalid deterministic rows (failed E*, boundary E*, complex/negative phi, ODE failure, NaN/Inf).

Stratifies remaining rows into low / medium / high `average_phi` (or `average_delta_mu`).

Uses `oversampleFactor` (default 1.25–1.5) so some SSA invalid rows do not empty the final plot.

### Stage 3 — SSA/KL

Gillespie SSA at E0 and E1 only for selected rows.

Default practical settings: `N_stoch=50`, `t_end_stoch=25`, `nBins=10`.

**Do not use `t_end_stoch=1500` as default** until SSA convergence is verified. Legacy value often hits `max_steps=1e6`.

Outputs: KL, symmetric KL, supplementary JS divergence, `valid_for_kl`, trajectory diagnostics.

### Stage 4 — plotting

Primary KL trend plots use **`valid_for_kl=true` only**. Invalid rows may appear as faint gray points.

## Example commands

### Smoke (~5–10 min)

```matlab
run('scripts/run_smoke_pipeline.m')
```

### First result (~30–90 min)

```matlab
run('scripts/run_first_result_pipeline.m')
```

### Manual stage control

```matlab
setup_repo_paths();
run_deterministic_screen('N_runs',500,'n_E_vals',100,'useParallel',true, ...
    'outputFile','results/first_result/results_deterministic_screen.csv','randomSeed',1);
select_parameter_sets('inputFile','results/first_result/results_deterministic_screen.csv', ...
    'N_select',50,'oversampleFactor',1.5,'randomSeed',2);
run_ssa_kl_on_selected('inputFile','results/first_result/selected_parameter_sets.csv', ...
    'N_stoch',50,'t_end_stoch',25,'nBins',10,'randomSeed',3);
plot_pipeline_results('outputDir','figures/first_result');
```

## Parameter ranges (unchanged from legacy)

- k1  = 10^(3 + 4·rand)
- k_1 = 10^(-1 + 4·rand)
- k2  = 10^(3 + 4·rand)
- k_2 = 10^(-1 + 4·rand)
- k3  = 10^(-1 + 4·rand)
- k_3 = 10^(1 + 4·rand)

## Input signal caveat

E0/E1 are currently defined by varying **total enzyme concentration**. If the intended input signal is substrate S or influx alpha, the scientific definition should be changed before final interpretation.

## What was not changed

- ODE equations and Gillespie reactions
- E* / E0 / E1 definitions
- KL histogram formula (linspace bins, PDF normalization, epsilon smoothing)
- Parameter sampling ranges

## Diagnostics

```matlab
setup_repo_paths();
check_metric_refactor_equivalence(5);
```

Verifies shared deterministic metrics match legacy separate ODE helpers.
