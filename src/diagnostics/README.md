# Diagnostics

Optional checks (not required for production pipeline runs).

## Metric equivalence

```matlab
setup_repo_paths();
check_metric_refactor_equivalence(5);
```

Confirms `compute_steady_state_metrics` matches legacy separate ODE helpers for phi, delta mu, and Jnet.

## Legacy timing diagnostic

The original timed diagnostic for the brute-force script is preserved in upstream development notes under `legacy/REVIEW_NOTES.md`. For SSA feasibility, use the first-result pipeline with `summary.md` validity counts before increasing `t_end_stoch` or `N_stoch`.
