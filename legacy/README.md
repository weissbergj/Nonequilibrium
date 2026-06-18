# Legacy code

Original upstream files preserved unchanged for reference.

| File | Description |
|------|-------------|
| `10k SAimulations` | Monolithic 10k parfor entry point |
| `Overallfinding3.m` | Valid MATLAB filename copy |
| `overall code` | Early ODE sweep prototype |
| `kl divergence`, `chemical potential`, `energy dissipiation` | Snippet files |
| `probability distribution`, `Total enzyme concentration` | Related snippets |
| `REVIEW_NOTES.md` | Refactor review notes from development |

**Do not use these as the main workflow.** Use `scripts/run_smoke_pipeline.m` or `scripts/run_first_result_pipeline.m` instead.

The staged pipeline in `src/` preserves the same model formulas but adds validity filtering, CSV outputs, and practical default SSA settings.
